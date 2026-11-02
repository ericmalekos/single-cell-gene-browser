"""GeneIndex resolution on small synthetic tables."""
from scbrowser.genes import CROSS_SPECIES_NOTE, GeneIndex
from scbrowser.ids import is_ensembl, split_ids, strip_ens


def test_strip_ens():
    assert strip_ens("ENSMUSG00000026104.16") == "ENSMUSG00000026104"
    assert strip_ens("ensg00000002586.20_PAR_Y") == "ENSG00000002586"
    assert strip_ens("Tex19.1") == "Tex19.1"          # dots in real symbols are kept
    assert is_ensembl(" ENSG00000141510 ") and not is_ensembl("H1-1")
    assert split_ids("ENSMUSG00000000001.5,ENSMUSG00000000003|") == ["ENSMUSG00000000001", "ENSMUSG00000000003"]


def mouse_index(**kw):
    symbols = ["Stat1", "Ptp4a1", "Ptp4a1", "Tex19.1", "ERCC-00002", "March1", "D17H6S56E-5"]
    ids = [["ENSMUSG00000026104"], ["ENSMUSG00000026064"], ["ENSMUSG00000099999"], [], [], ["ENSMUSG00000036469"],
           ["ENSMUSG00000120992"]]
    return GeneIndex(symbols, ids, species="mouse", hidden=[False, False, False, False, True, False, False],
                     current=[None, None, None, None, None, "Marchf1", None],
                     aliases={"synir": ["D17H6S56E-5"], "isgf-3": ["Stat1"]}, **kw)


def test_symbol_and_case():
    idx = mouse_index()
    r = idx.resolve("Stat1")
    assert r.rows == [0] and r.how == "symbol" and not r.note
    r = idx.resolve("stat1")
    assert r.rows == [0] and "case-insensitive" in r.how and not r.note   # lower case is not cross-species


def test_human_style_query_on_mouse_gets_note():
    r = mouse_index().resolve("STAT1")
    assert r.rows == [0] and r.note == CROSS_SPECIES_NOTE


def test_ensembl_versioned_and_not():
    idx = mouse_index()
    assert idx.resolve("ENSMUSG00000026104").rows == [0]
    assert idx.resolve("ensmusg00000026104.16").rows == [0]
    assert idx.resolve("Tex19.1").rows == [3]


def test_duplicates_and_outdated_symbols():
    idx = mouse_index()
    assert idx.resolve("Ptp4a1").rows == [1, 2]
    r = idx.resolve("Marchf1")                      # dataset lists the gene under its old symbol
    assert r.rows == [5] and "March1" in r.how


def test_aliases():
    idx = mouse_index()
    assert idx.resolve("Synir").rows == [6]
    assert idx.resolve("ISGF-3").rows == [0]


def test_hidden_rows_searchable_but_not_completed():
    idx = mouse_index()
    assert idx.resolve("ERCC-00002").rows == [4]
    assert "ERCC-00002" not in idx.completions() and "Stat1" in idx.completions()


def test_not_found_suggestions():
    r = mouse_index().resolve("Stat")
    assert not r.found and r.suggestions[0] == "Stat1"
    r = mouse_index().resolve("Ptp4a2")             # close spelling
    assert "Ptp4a1" in r.suggestions


def test_make_unique_variants_and_foreign_ids():
    symbols = ["MATR3", "MATR3-1", "STAT1", "H1-1"]
    idx = GeneIndex(symbols, [["ENSG00000015479", "ENSG00000280987"], ["ENSG00000015479", "ENSG00000280987"],
                              ["ENSG00000115415"], ["ENSG00000124610"]],
                    species="human", variant_of=[None, "MATR3", None, float("nan")],
                    foreign_ids=lambda gid: [("Stat1", "mouse")] if gid == "ENSMUSG00000026104" else [])
    assert idx.resolve("MATR3").rows == [0, 1]
    assert idx.resolve("MATR3-1").rows == [1, 0]
    assert idx.resolve("H1-1").rows == [3]
    r = idx.resolve("ENSMUSG00000026104.16")
    assert r.rows == [2] and r.note == CROSS_SPECIES_NOTE
    assert not idx.resolve("ENSMUSG00000000001").found
