"""
VCMC v3 — Decoder
from scratch, by ananya

Reverses the VCMC cipher. Given a SMILES string or a docked PDB/PDBQT file,
outputs the original word that was encoded.

Usage:
    python3 vcmc_decode.py "SMILES_STRING"
    python3 vcmc_decode.py cant_vcmc_docked.pdbqt
    python3 vcmc_decode.py --folder ./docked/

The decoder reads atom-by-atom, identifies substituent + chirality,
and maps back to the original letter using the reverse LETTER_MAP.
VIBGYOR positional tag (Rule 05) is used to reconstruct sentence word order.
"""

import sys
import os
import glob
from rdkit import Chem

# ── Forward map (same as encoder) ─────────────────────────────
LETTER_MAP = {
    'A': 'N',        'E': 'O',        'I': 'N(C)',
    'O': 'S',        'U': 'C(=O)',
    'B': '[C@@H]',   'C': '[C@H]',
    'D': '[C@@H](F)', 'F': '[C@H](F)',
    'G': '[C@@H](Cl)','H': '[C@H](Cl)',
    'J': '[C@@H](Br)','K': '[C@H](Br)',
    'L': '[C@@H](O)', 'M': '[C@H](O)',
    'V': '[C@@H](S)', 'W': '[C@H](S)',
    'N': '[C@@H](N)', 'P': '[C@H](N)',
    'Q': '[C@@H](C)', 'R': '[C@H](C)',
    'S': '[C@@H](CC)','T': '[C@H](CC)',
    'X': '[C@@H](C#N)','Y': '[C@H](C#N)',
    'Z': '[C@@H](I)',
}

# ── Reverse map: SMILES tag → letter ──────────────────────────
REVERSE_MAP = {v: k for k, v in LETTER_MAP.items()}

# ── VIBGYOR sentence order colours ────────────────────────────
VIBGYOR = ['Violet','Indigo','Blue','Green','Yellow','Orange','Red']

def smiles_to_word(smiles):
    """Decode a SMILES string back to the original word."""
    smiles = smiles.strip()

    # Canonicalise for consistent parsing
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  ✗ Could not parse SMILES: {smiles}")
        return None

    # Re-generate SMILES with stereo info preserved
    # We parse the raw SMILES token by token instead,
    # since RDKit canonicalisation can reorder atoms

    # ── Token parser ──────────────────────────────────────────
    # Strips ring closure digits and parses atom tags in order
    import re

    # Remove ring closure digits (e.g. [C@@H](CC)1...1 → remove the 1s)
    raw = re.sub(r'(?<=[A-Za-z\]])(\d)', '', smiles)

    # Match atom tokens in order: [C@@H](XX), N(C), C(=O), N, O, S
    pattern = re.compile(
        r'(\[C@@H\]\(C#N\)|\[C@H\]\(C#N\)'   # nitrile
        r'|\[C@@H\]\(CC\)|\[C@H\]\(CC\)'      # ethyl
        r'|\[C@@H\]\(Cl\)|\[C@H\]\(Cl\)'      # chloro
        r'|\[C@@H\]\(Br\)|\[C@H\]\(Br\)'      # bromo
        r'|\[C@@H\]\(F\)|\[C@H\]\(F\)'        # fluoro
        r'|\[C@@H\]\(O\)|\[C@H\]\(O\)'        # hydroxyl
        r'|\[C@@H\]\(S\)|\[C@H\]\(S\)'        # thiol
        r'|\[C@@H\]\(N\)|\[C@H\]\(N\)'        # amino
        r'|\[C@@H\]\(C\)|\[C@H\]\(C\)'        # methyl
        r'|\[C@@H\]\(I\)|\[C@H\]\(I\)'        # iodo
        r'|\[C@@H\]|\[C@H\]'                  # bare chiral C
        r'|N\(C\)|C\(=O\)'                    # vowels with groups
        r'|[NOSC])'                            # plain vowel atoms
    )

    tokens = pattern.findall(raw)

    if not tokens:
        print(f"  ✗ No recognisable VCMC atoms found in: {smiles}")
        return None

    word = ''
    print(f"\n  {'ATOM TAG':<18} {'LETTER':<8} {'ROLE'}")
    print(f"  {'-'*50}")

    for token in tokens:
        letter = REVERSE_MAP.get(token)
        if letter:
            role = 'vowel → heteroatom' if letter in 'AEIOU' else 'consonant → chiral C'
            print(f"  {token:<18} {letter:<8} {role}")
            word += letter
        else:
            print(f"  {token:<18} {'?':<8} unrecognised — skipped")

    return word


def decode_pdb_or_pdbqt(filepath):
    """Extract SMILES from a PDB/PDBQT file and decode."""
    # Try reading with RDKit
    if filepath.endswith('.pdbqt'):
        # Convert PDBQT → PDB on the fly by stripping extra columns
        lines = []
        with open(filepath) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM', 'END', 'MODEL', 'ENDMDL')):
                    lines.append(line[:66])  # standard PDB width
        from rdkit.Chem import MolFromPDBBlock
        pdb_block = ''.join(lines)
        mol = MolFromPDBBlock(pdb_block, removeHs=False)
    else:
        mol = Chem.MolFromPDBFile(filepath, removeHs=False)

    if mol is None:
        print(f"  ✗ Could not read: {filepath}")
        return None

    # Get SMILES with stereo
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    print(f"  Extracted SMILES: {smiles}")
    return smiles_to_word(smiles)


def decode_folder(folder):
    """Decode all docked PDBQT files in a folder, ordered by VIBGYOR position."""
    files = sorted(glob.glob(os.path.join(folder, '*_docked.pdbqt')))
    if not files:
        files = sorted(glob.glob(os.path.join(folder, '*_vcmc.pdb')))
    if not files:
        print(f"No docked PDBQT or PDB files found in {folder}")
        return

    print(f"\n{'='*60}")
    print(f"  VCMC v3 · Sentence Decoder")
    print(f"  {len(files)} molecule(s) found")
    print(f"{'='*60}")

    sentence = []
    for i, fpath in enumerate(files):
        basename = os.path.basename(fpath)
        colour = VIBGYOR[i] if i < len(VIBGYOR) else f'pos-{i+1}'
        print(f"\n[{colour}] File: {basename}")
        print(f"  {'─'*50}")
        word = decode_pdb_or_pdbqt(fpath)
        if word:
            sentence.append((colour, word))
            print(f"\n  ✓ Decoded word: {word}")

    print(f"\n{'='*60}")
    print(f"  DECODED SENTENCE:")
    print(f"  {' '.join(w for _, w in sentence)}")
    print(f"\n  Word order (Rule 05 — VIBGYOR):")
    for colour, word in sentence:
        print(f"    {colour:<10} → {word}")
    print(f"{'='*60}\n")


def decode_smiles_string(smiles):
    """Decode a raw SMILES string."""
    print(f"\n{'='*60}")
    print(f"  VCMC v3 · SMILES Decoder")
    print(f"  Input: {smiles}")
    print(f"{'='*60}")
    word = smiles_to_word(smiles)
    if word:
        print(f"\n  ✓ Decoded word: {word}")
    print(f"{'='*60}\n")
    return word


# ── Entry point ───────────────────────────────────────────────
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    arg = sys.argv[1]

    if arg == '--folder' and len(sys.argv) > 2:
        decode_folder(sys.argv[2])
    elif os.path.isdir(arg):
        decode_folder(arg)
    elif os.path.isfile(arg):
        print(f"\n{'='*60}")
        print(f"  VCMC v3 · File Decoder")
        print(f"  File: {arg}")
        print(f"{'='*60}")
        word = decode_pdb_or_pdbqt(arg)
        if word:
            print(f"\n  ✓ Decoded word: {word}")
        print(f"{'='*60}\n")
    else:
        # Treat as SMILES string
        decode_smiles_string(arg)
