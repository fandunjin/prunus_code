"""EDTA class and REXdb subfamily labels used in exclusive coverage."""

import re


def normalize_class(token: str) -> str | None:
    """Normalize an EDTA GFF classification token."""
    t = token.lower().replace("_", "-")
    if "gypsy" in t:
        return "LTR/Gypsy"
    if "copia" in t:
        return "LTR/Copia"
    if t.startswith("ltr/"):
        return "Other LTR"
    if "helitron" in t:
        return "Helitron"
    if t.startswith("tir/") or t.startswith("dna/"):
        return "TIR/DNA"
    if t.startswith("line"):
        return "LINE"
    if t.startswith("sine") or "trna" in t:
        return "SINE"
    if t in {"snrna/unknown", "rrna/unknown", "repeat-fragment", "repeat_fragment"}:
        return "Other TE"
    if t.startswith("nonltr/") or t.startswith("non-tir/") or t.startswith("nonltr"):
        return "Other TE"
    if t.startswith("unknown/") or t == "unknown":
        return "Other TE"
    if "rdna" in t or t.startswith("45s") or t.startswith("rrna"):
        return None
    return "Other TE" if "/" in t else None


def dna_subfamily(token: str) -> str:
    t = token.lower()
    if "helitron" in t:
        return "Helitron"
    if any(k in t for k in ("dtm", "mule", "mutator")):
        return "Mutator"
    if any(k in t for k in ("dta", "hat")):
        return "hAT"
    if any(k in t for k in ("dtc", "cmc", "enspm", "cacta")):
        return "CACTA"
    if any(k in t for k in ("dth", "pif", "harbinger", "tourist")):
        return "PIF/Harbinger"
    if any(k in t for k in ("dtt", "tc1", "mariner")):
        return "Tc1/Mariner"
    return "other TIR/DNA"


def line_subfamily(token: str) -> str:
    t = token.lower()
    if t.startswith("line/l1"):
        return "L1"
    if "rte" in t:
        return "RTE"
    return "other LINE"


def rexdb_lineage(family_cls: str, info: tuple[str, ...]) -> str:
    """info = (class, order, family, sublineage, clade) from REXdb v3.0.

    For Ty3/gypsy rows the 6th column carry chromovirus/non-chromovirus and the
    clade (OTA, CRM, Tekay, ...). For Ty1/copia rows the lineage name
    (Ale/SIRE/Ivana/Tork/...) sits in the 4th column (sublineage)."""
    _cls, _order, fam, col4, col5 = info
    if family_cls == "LTR/Gypsy":
        clade = col5 or col4
        if clade == "OTA":
            return "OTA"
        if clade == "CRM":
            return "CRM"
        if clade == "Tekay":
            return "Tekay"
        if col4 == "chromovirus" and clade not in (
            "",
            "chromo-unclass",
            "chromo-outgroup",
        ):
            return "other chromovirus"
        return "other Gypsy"
    if family_cls == "LTR/Copia":
        lineage = col4 or col5
        if lineage in ("Ale", "SIRE", "Ivana", "Tork"):
            return lineage
        return "other Copia"
    if family_cls == "Other LTR":
        if fam == "Ty3/gypsy":
            return "Gypsy-like"
        if fam == "Ty1/copia":
            return "Copia-like"
        return "unclassified"
    return ""


def subfamily_of(cls: str, token: str, fam: str, lin: dict) -> str:
    if cls in ("LTR/Gypsy", "LTR/Copia", "Other LTR"):
        info = lin.get(fam)
        if info is None and "_" in fam:
            base = re.sub(r"_(INT|LTR)$", "", fam)
            info = lin.get(base)
        if info:
            lab = rexdb_lineage(cls, info)
            if lab:
                return lab
        return (
            "other Gypsy"
            if cls == "LTR/Gypsy"
            else ("other Copia" if cls == "LTR/Copia" else "unclassified")
        )
    if cls == "TIR/DNA":
        return dna_subfamily(token)
    if cls == "LINE":
        return line_subfamily(token)
    return cls
