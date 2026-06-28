"""
VCMC v3 — Decoder
from scratch, by ananya

Reads SMILES directly from *_vcmc.pdb files (stereo preserved).
Reconstructs original words and sentence in VIBGYOR order.

Usage:
    python3 vcmc_decode.py                        # decode all *_vcmc.pdb in current folder
    python3 vcmc_decode.py --folder ./some/path/  # decode from specific folder
    python3 vcmc_decode.py "[C@H]N[C@@H](N)[C@H](CC)"  # decode single SMILES
"""

import sys, os, glob, re
from rdkit import Chem

LETTER_MAP = {
    'A':'N','E':'O','I':'N(C)','O':'S','U':'C(=O)',
    'B':'[C@@H]','C':'[C@H]',
    'D':'[C@@H](F)','F':'[C@H](F)',
    'G':'[C@@H](Cl)','H':'[C@H](Cl)',
    'J':'[C@@H](Br)','K':'[C@H](Br)',
    'L':'[C@@H](O)','M':'[C@H](O)',
    'V':'[C@@H](S)','W':'[C@H](S)',
    'N':'[C@@H](N)','P':'[C@H](N)',
    'Q':'[C@@H](C)','R':'[C@H](C)',
    'S':'[C@@H](CC)','T':'[C@H](CC)',
    'X':'[C@@H](C#N)','Y':'[C@H](C#N)',
    'Z':'[C@@H](I)',
}
REVERSE_MAP = {v: k for k, v in LETTER_MAP.items()}
VIBGYOR = ['Violet','Indigo','Blue','Green','Yellow','Orange','Red']

TOKEN_RE = re.compile(
    r'\[C@@H\]\(C#N\)|\[C@H\]\(C#N\)'
    r'|\[C@@H\]\(CC\)|\[C@H\]\(CC\)'
    r'|\[C@@H\]\(Cl\)|\[C@H\]\(Cl\)'
    r'|\[C@@H\]\(Br\)|\[C@H\]\(Br\)'
    r'|\[C@@H\]\(F\)|\[C@H\]\(F\)'
    r'|\[C@@H\]\(O\)|\[C@H\]\(O\)'
    r'|\[C@@H\]\(S\)|\[C@H\]\(S\)'
    r'|\[C@@H\]\(N\)|\[C@H\]\(N\)'
    r'|\[C@@H\]\(C\)|\[C@H\]\(C\)'
    r'|\[C@@H\]\(I\)|\[C@H\]\(I\)'
    r'|\[C@@H\]|\[C@H\]'
    r'|N\(C\)|C\(=O\)'
    r'|(?<![A-Za-z])[NOSA](?![A-Za-z])'
)

def smiles_to_word(smiles):
    raw = re.sub(r'(?<=[A-Za-z\]])(\d)', '', smiles.strip())
    tokens = TOKEN_RE.findall(raw)
    if not tokens:
        return None
    word = ''
    print(f"\n  {'ATOM TAG':<20} {'LETTER':<8} ROLE")
    print(f"  {'-'*52}")
    for token in tokens:
        letter = REVERSE_MAP.get(token)
        if letter:
            role = 'vowel → heteroatom' if letter in 'AEIOU' else 'consonant → chiral C'
            print(f"  {token:<20} {letter:<8} {role}")
            word += letter
        else:
            print(f"  {token:<20} {'?':<8} unrecognised — skipped")
    return word

def get_smiles_from_pdb(filepath):
    """Extract SMILES with stereo from a *_vcmc.pdb file."""
    mol = Chem.MolFromPDBFile(filepath, removeHs=True, sanitize=True)
    if mol is None:
        return None
    # isomericSmiles=True preserves [C@@H] / [C@H] stereochemistry
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False)
    return smiles

def decode_folder(folder):
    # Always use *_vcmc.pdb — stereo is fully preserved there
    files = sorted(glob.glob(os.path.join(folder, '*_vcmc.pdb')))
    if not files:
        print(f"No *_vcmc.pdb files found in {folder}")
        print("Make sure you're pointing at the folder with the original encoded PDB files.")
        return

    print(f"\n{'='*60}")
    print(f"  VCMC v3 · Sentence Decoder")
    print(f"  {len(files)} word(s) found")
    print(f"{'='*60}")

    sentence = []
    for i, fpath in enumerate(files):
        basename = os.path.basename(fpath)
        colour = VIBGYOR[i] if i < len(VIBGYOR) else f'pos-{i+1}'
        print(f"\n[{colour}] {basename}")
        print(f"  {'─'*50}")
        smiles = get_smiles_from_pdb(fpath)
        if smiles:
            print(f"  SMILES: {smiles}")
            word = smiles_to_word(smiles)
            if word:
                sentence.append((colour, word))
                print(f"\n  ✓ Decoded: {word}")
        else:
            print(f"  ✗ Could not read {basename}")

    print(f"\n{'='*60}")
    print(f"  DECODED SENTENCE: {' '.join(w for _, w in sentence)}")
    print(f"\n  VIBGYOR word order:")
    for colour, word in sentence:
        print(f"    {colour:<10} → {word}")
    print(f"{'='*60}\n")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        decode_folder('.')
    elif sys.argv[1] == '--folder' and len(sys.argv) > 2:
        decode_folder(sys.argv[2])
    elif os.path.isdir(sys.argv[1]):
        decode_folder(sys.argv[1])
    else:
        # treat as SMILES
        smiles = sys.argv[1]
        print(f"\n{'='*60}\n  VCMC v3 · SMILES Decoder\n  Input: {smiles}\n{'='*60}")
        word = smiles_to_word(smiles)
        if word:
            print(f"\n  ✓ Decoded: {word}")
        print(f"{'='*60}\n")
