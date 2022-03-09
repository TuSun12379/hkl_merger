
STBF = 12    #  Seitz Matrix Translation          Base Factor */
CRBF = 12    # Change of Basis Matrix Rotation    Base Factor */
CTBF = 72    # Change of Basis Matrix Translation Base Factor */
SRBF = 1

PG_Names = ["Unknown", "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4",
    "4/m", "422", "4mm", "-42m", "-4m2", "4/mmm", "3", "-3", "321", "312", "32",
    "3m1", "31m", "3m", "-3m1", "-31m", "-3m", "6", "-6", "6/m", "622", "6mm", "-6m2",
    "-62m", "6/mmm", "23", "m-3", "432", "-43m", "m-3m"]

XS_Name = ["Unknown", "Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal", "Trigonal", "Hexagonal", "Cubic"]

EI_Name = ["Unknown", "Enantiomorphic", "Obverse", "Reverse"]

SchoenfliesSymbols = [None, "C1^1", "Ci^1", "C2^1", "C2^2", "C2^3", "Cs^1", "Cs^2", "Cs^3", "Cs^4", "C2h^1", "C2h^2",
    "C2h^3", "C2h^4", "C2h^5", "C2h^6", "D2^1", "D2^2", "D2^3", "D2^4", "D2^5", "D2^6", "D2^7", "D2^8", "D2^9", "C2v^1",
    "C2v^2", "C2v^3", "C2v^4", "C2v^5", "C2v^6", "C2v^7", "C2v^8", "C2v^9", "C2v^10", "C2v^11", "C2v^12", "C2v^13", "C2v^14",
    "C2v^15", "C2v^16", "C2v^17", "C2v^18", "C2v^19", "C2v^20", "C2v^21", "C2v^22", "D2h^1", "D2h^2", "D2h^3", "D2h^4", "D2h^5",
    "D2h^6", "D2h^7", "D2h^8", "D2h^9", "D2h^10", "D2h^11", "D2h^12", "D2h^13", "D2h^14", "D2h^15", "D2h^16", "D2h^17", "D2h^18",
    "D2h^19", "D2h^20", "D2h^21", "D2h^22", "D2h^23", "D2h^24", "D2h^25", "D2h^26", "D2h^27", "D2h^28", "C4^1", "C4^2", "C4^3",
    "C4^4", "C4^5", "C4^6", "S4^1", "S4^2", "C4h^1", "C4h^2", "C4h^3", "C4h^4", "C4h^5", "C4h^6", "D4^1", "D4^2", "D4^3", "D4^4",
    "D4^5", "D4^6", "D4^7", "D4^8", "D4^9", "D4^10", "C4v^1", "C4v^2", "C4v^3", "C4v^4", "C4v^5", "C4v^6", "C4v^7", "C4v^8", "C4v^9",
    "C4v^10", "C4v^11", "C4v^12", "D2d^1", "D2d^2", "D2d^3", "D2d^4", "D2d^5", "D2d^6", "D2d^7", "D2d^8", "D2d^9", "D2d^10", "D2d^11",
    "D2d^12", "D4h^1", "D4h^2", "D4h^3", "D4h^4", "D4h^5", "D4h^6", "D4h^7", "D4h^8", "D4h^9", "D4h^10", "D4h^11", "D4h^12", "D4h^13",
    "D4h^14", "D4h^15", "D4h^16", "D4h^17", "D4h^18", "D4h^19", "D4h^20", "C3^1", "C3^2", "C3^3", "C3^4", "C3i^1", "C3i^2", "D3^1",
    "D3^2", "D3^3", "D3^4", "D3^5", "D3^6", "D3^7", "C3v^1", "C3v^2", "C3v^3", "C3v^4", "C3v^5", "C3v^6", "D3d^1", "D3d^2", "D3d^3",
    "D3d^4", "D3d^5", "D3d^6", "C6^1", "C6^2", "C6^3", "C6^4", "C6^5", "C6^6", "C3h^1", "C6h^1", "C6h^2", "D6^1", "D6^2", "D6^3",
    "D6^4", "D6^5", "D6^6", "C6v^1", "C6v^2", "C6v^3", "C6v^4", "D3h^1", "D3h^2", "D3h^3", "D3h^4", "D6h^1", "D6h^2", "D6h^3",
    "D6h^4", "T^1", "T^2", "T^3", "T^4", "T^5", "Th^1", "Th^2", "Th^3", "Th^4", "Th^5", "Th^6", "Th^7", "O^1", "O^2", "O^3",
    "O^4", "O^5", "O^6", "O^7", "O^8", "Td^1", "Td^2", "Td^3", "Td^4", "Td^5", "Td^6", "Oh^1", "Oh^2", "Oh^3", "Oh^4", "Oh^5",
    "Oh^6", "Oh^7", "Oh^8", "Oh^9", "Oh^10"]
