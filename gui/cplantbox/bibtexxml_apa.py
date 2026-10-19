"""Convert bibliography XML records to APA-like citation strings.

Supports both:
- BibTeXML style records (<entry id="..."> ...)
- MODS XML records produced by bibutils/bib2xml (<modsCollection> / <mods>)

Typical workflow:
1. Convert .bib to XML, e.g. bib2xml myrefs.bib > myrefs.xml
2. Run this module on the XML file to extract and format citations.
"""

import re
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

# from __future__ import annotations


def _local_name(tag: str) -> str:
    """Return XML local tag name without namespace."""
    return tag.split("}", 1)[-1] if "}" in tag else tag


def _normalize_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _collect_text(element: ET.Element) -> str:
    """Collect all text recursively from an element."""
    return _normalize_spaces("".join(element.itertext()))


def _extract_fields(entry_type_node: ET.Element) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    authors: List[str] = []

    for child in list(entry_type_node):
        key = _local_name(child.tag).lower()
        value = _collect_text(child)
        if not value:
            continue

        if key == "author":
            authors.append(value)
            continue

        if key in fields and fields[key]:
            fields[key] = f"{fields[key]} and {value}"
        else:
            fields[key] = value

    if authors:
        fields["author"] = " and ".join(authors)

    return fields


def _children_by_local(parent: ET.Element, name: str) -> List[ET.Element]:
    """Return direct children whose local tag name matches ``name``."""
    lname = name.lower()
    return [child for child in list(parent) if _local_name(child.tag).lower() == lname]


def _first_child_by_local(parent: ET.Element, name: str) -> Optional[ET.Element]:
    """Return first direct child by local name, or None."""
    for child in list(parent):
        if _local_name(child.tag).lower() == name.lower():
            return child
    return None


def _extract_mods_entries(root: ET.Element) -> List[Dict[str, str]]:
    """Extract citation-like entries from MODS XML (e.g. bibutils bib2xml output)."""
    entries: List[Dict[str, str]] = []

    for node in root.iter():
        if _local_name(node.tag).lower() != "mods":
            continue

        fields: Dict[str, str] = {}

        # title + optional subtitle
        title_info = _first_child_by_local(node, "titleInfo")
        if title_info is not None:
            title_node = _first_child_by_local(title_info, "title")
            subtitle_node = _first_child_by_local(title_info, "subTitle")
            title = _collect_text(title_node) if title_node is not None else ""
            subtitle = _collect_text(subtitle_node) if subtitle_node is not None else ""
            if title and subtitle:
                fields["title"] = f"{title}: {subtitle}"
            elif title:
                fields["title"] = title

        # authors
        authors: List[str] = []
        for name_node in _children_by_local(node, "name"):
            role = ""
            role_node = _first_child_by_local(name_node, "role")
            if role_node is not None:
                role_term = _first_child_by_local(role_node, "roleTerm")
                role = _collect_text(role_term).lower() if role_term is not None else ""

            if role and role != "author":
                continue

            family = ""
            givens: List[str] = []
            for part in _children_by_local(name_node, "namePart"):
                ptype = (part.attrib.get("type") or "").lower()
                value = _collect_text(part)
                if not value:
                    continue
                if ptype == "family":
                    family = value
                elif ptype == "given":
                    givens.append(value)

            if family and givens:
                authors.append(f"{family}, {' '.join(givens)}")
            elif family:
                authors.append(family)
            elif givens:
                authors.append(" ".join(givens))

        if authors:
            fields["author"] = " and ".join(authors)

        # publication year
        year = ""
        origin = _first_child_by_local(node, "originInfo")
        if origin is not None:
            date_issued = _first_child_by_local(origin, "dateIssued")
            if date_issued is not None:
                year = _collect_text(date_issued)
        if not year:
            part = _first_child_by_local(node, "part")
            if part is not None:
                date = _first_child_by_local(part, "date")
                if date is not None:
                    year = _collect_text(date)
        if year:
            fields["year"] = year

        # journal + publisher (host item)
        related_host = None
        for rel in _children_by_local(node, "relatedItem"):
            if (rel.attrib.get("type") or "").lower() == "host":
                related_host = rel
                break

        if related_host is not None:
            host_title_info = _first_child_by_local(related_host, "titleInfo")
            if host_title_info is not None:
                host_title = _first_child_by_local(host_title_info, "title")
                if host_title is not None:
                    journal = _collect_text(host_title)
                    if journal:
                        fields["journal"] = journal

            host_origin = _first_child_by_local(related_host, "originInfo")
            if host_origin is not None:
                publisher = _first_child_by_local(host_origin, "publisher")
                if publisher is not None:
                    pub_val = _collect_text(publisher)
                    if pub_val:
                        fields["publisher"] = pub_val

        # pages (bibutils often maps article id/page-like token to this)
        part = _first_child_by_local(node, "part")
        if part is not None:
            for detail in _children_by_local(part, "detail"):
                if (detail.attrib.get("type") or "").lower() == "page":
                    number = _first_child_by_local(detail, "number")
                    if number is not None:
                        pages = _collect_text(number)
                        if pages:
                            fields["pages"] = pages
                    break

        # citation key
        entry_id = node.attrib.get("ID", "")
        if not entry_id:
            for ident in _children_by_local(node, "identifier"):
                if (ident.attrib.get("type") or "").lower() == "citekey":
                    entry_id = _collect_text(ident)
                    if entry_id:
                        break

        entries.append({"id": entry_id, "type": "article", "fields": fields})

    return entries


def _initials(given_names: str) -> str:
    tokens = [t for t in re.split(r"[\s-]+", given_names.strip()) if t]
    return " ".join(f"{t[0].upper()}." for t in tokens)


def _format_one_author(raw_author: str) -> str:
    author = _normalize_spaces(raw_author)
    if not author:
        return ""

    if "," in author:
        last, given = [p.strip() for p in author.split(",", 1)]
    else:
        parts = author.split(" ")
        if len(parts) == 1:
            return parts[0]
        last = parts[-1]
        given = " ".join(parts[:-1])

    given_initials = _initials(given)
    return f"{last}, {given_initials}" if given_initials else last


def _format_authors_apa(author_field: str) -> str:
    people = [_format_one_author(a) for a in author_field.split(" and ")]
    people = [p for p in people if p]

    if not people:
        return ""
    if len(people) == 1:
        return people[0]
    if len(people) == 2:
        return f"{people[0]}, & {people[1]}"
    return ", ".join(people[:-1]) + f", & {people[-1]}"


def _clean_title(title: str) -> str:
    title = _normalize_spaces(title)
    if not title:
        return ""
    return title[:-1] if title.endswith(".") else title


def _format_doi(doi: str) -> str:
    doi = _normalize_spaces(doi)
    if not doi:
        return ""
    if doi.startswith("http://") or doi.startswith("https://"):
        return doi
    return f"https://doi.org/{doi}"


def _format_doi_markdown(doi: str) -> str:
    """Return DOI as markdown hyperlink text."""
    doi_url = _format_doi(doi)
    if not doi_url:
        return ""
    return f"[DOI]({doi_url})"


def bibtex_item_to_bibtexxml(bibtex_entry: str) -> str:
    """Convert one BibTeX item string (e.g., @article{key, ...}) to a BibTeXML <entry> string."""

    text = bibtex_entry.strip()
    match = re.match(r"@\s*(\w+)\s*\{\s*([^,\s]+)\s*,(.*)\}\s*$", text, re.DOTALL)
    if not match:
        raise ValueError("Input is not a valid single BibTeX entry")

    entry_type, entry_id, body = match.groups()
    fields: Dict[str, str] = {}

    i = 0
    n = len(body)
    while i < n:
        while i < n and body[i] in " \t\r\n,":
            i += 1
        if i >= n:
            break

        key_start = i
        while i < n and (body[i].isalnum() or body[i] in "_-:"):
            i += 1
        key = body[key_start:i].strip().lower()
        if not key:
            break

        while i < n and body[i].isspace():
            i += 1
        if i >= n or body[i] != "=":
            raise ValueError(f"Expected '=' after field name '{key}'")
        i += 1

        while i < n and body[i].isspace():
            i += 1
        if i >= n:
            break

        if body[i] == "{":
            depth = 1
            i += 1
            value_start = i
            while i < n and depth > 0:
                if body[i] == "{":
                    depth += 1
                elif body[i] == "}":
                    depth -= 1
                i += 1
            if depth != 0:
                raise ValueError(f"Unbalanced braces in field '{key}'")
            value = body[value_start : i - 1]
        elif body[i] == '"':
            i += 1
            value_start = i
            escaped = False
            while i < n:
                ch = body[i]
                if ch == '"' and not escaped:
                    break
                escaped = (ch == "\\") and not escaped
                if ch != "\\":
                    escaped = False
                i += 1
            if i >= n or body[i] != '"':
                raise ValueError(f"Unbalanced quotes in field '{key}'")
            value = body[value_start:i]
            i += 1
        else:
            value_start = i
            while i < n and body[i] != ",":
                i += 1
            value = body[value_start:i]

        fields[key] = _normalize_spaces(value)

        while i < n and body[i] in " \t\r\n,":
            i += 1

    entry_elem = ET.Element("entry", {"id": entry_id})
    type_elem = ET.SubElement(entry_elem, entry_type.lower())

    for key, value in fields.items():
        if key == "author":
            for author in [a.strip() for a in value.split(" and ") if a.strip()]:
                author_elem = ET.SubElement(type_elem, "author")
                author_elem.text = author
        else:
            field_elem = ET.SubElement(type_elem, key)
            field_elem.text = value

    return ET.tostring(entry_elem, encoding="unicode")


def _format_apa(entry_type: str, fields: Dict[str, str]) -> str:
    authors = _format_authors_apa(fields.get("author", ""))
    year = fields.get("year", "n.d.")
    title = _clean_title(fields.get("title", "Untitled"))

    prefix = ""
    if authors:
        prefix = f"{authors} "
    prefix += f"({year}). {title}."

    entry_type = entry_type.lower()

    if entry_type == "article":
        journal = fields.get("journal", "")
        volume = fields.get("volume", "")
        number = fields.get("number", "")
        pages = fields.get("pages", "")
        doi_md = _format_doi_markdown(fields.get("doi", ""))

        details = ""
        if journal:
            details += f" *{journal}*"
        if volume:
            details += f", {volume}"
        if number:
            details += f"({number})"
        if pages:
            details += f", {pages}"
        if details:
            details += "."
        if doi_md:
            details += f" {doi_md}"
        return (prefix + details).strip()

    if entry_type == "book":
        publisher = fields.get("publisher", "")
        return f"{prefix} {publisher}.".strip() if publisher else prefix

    if entry_type == "inproceedings":
        booktitle = fields.get("booktitle", "")
        pages = fields.get("pages", "")
        publisher = fields.get("publisher", "")

        details = ""
        if booktitle:
            details += f" In {booktitle}"
        if pages:
            details += f" (pp. {pages})"
        if details:
            details += "."
        if publisher:
            details += f" {publisher}."
        return (prefix + details).strip()

    # Fallback for unsupported entry types.
    extras = []
    for key in ("journal", "booktitle", "publisher", "pages", "doi"):
        value = fields.get(key)
        if value:
            if key == "doi":
                extras.append(f"{key}={_format_doi_markdown(value)}")
            else:
                extras.append(f"{key}={value}")
    return (prefix + (" " + "; ".join(extras) if extras else "")).strip()


def extract_bibtexxml_entries(xml_file: str) -> List[Dict[str, str]]:
    """
    Extract BibTeXML entries from an XML file that may contain unrelated tags.

    Returns one dictionary per citation with fields:
    - id: citation key
    - type: bib entry type (article, book, inproceedings, ...)
    - fields: raw field dictionary
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    entries: List[Dict[str, str]] = []

    for node in root.iter():
        if _local_name(node.tag).lower() != "entry":
            continue

        entry_id = node.attrib.get("id", "")
        entry_type_node: Optional[ET.Element] = None
        for child in list(node):
            entry_type_node = child
            break

        if entry_type_node is None:
            continue

        entry_type = _local_name(entry_type_node.tag)
        fields = _extract_fields(entry_type_node)
        entries.append({"id": entry_id, "type": entry_type, "fields": fields})

    if entries:
        return entries

    # Fallback for bibutils MODS output (root is typically <modsCollection>).
    return _extract_mods_entries(root)


def bibtexxml_to_apa(xml_file: str, entry_id: Optional[str] = None) -> List[str]:
    """
    Convert BibTeXML entries from ``xml_file`` to APA-like citation strings.

    Parameters
    ----------
    xml_file
        Path to an XML file containing BibTeXML ``entry`` tags.
    entry_id
        Optional BibTeXML citation key to select a single entry.

    Returns
    -------
    List[str]
        One APA-like string per matched entry.
    """
    entries = extract_bibtexxml_entries(xml_file)

    if entry_id is not None:
        entries = [e for e in entries if e["id"] == entry_id]

    return [_format_apa(e["type"], e["fields"]) for e in entries]


def bibtexxml_apa_from_string(xml_str: str, entry_id: Optional[str] = None) -> List[str]:
    """
    Like bibtexxml_to_apa but accepts an XML string instead of a file path.

    Returns an empty list if the XML contains no recognisable citation entries.
    """
    try:
        root = ET.fromstring(xml_str)
    except ET.ParseError:
        return []

    entries: List[Dict[str, str]] = []
    for node in root.iter():
        if _local_name(node.tag).lower() != "entry":
            continue
        entry_id_ = node.attrib.get("id", "")
        entry_type_node: Optional[ET.Element] = None
        for child in list(node):
            entry_type_node = child
            break
        if entry_type_node is None:
            continue
        entry_type = _local_name(entry_type_node.tag)
        fields = _extract_fields(entry_type_node)
        entries.append({"id": entry_id_, "type": entry_type, "fields": fields})

    if not entries:
        entries = _extract_mods_entries(root)

    if entry_id is not None:
        entries = [e for e in entries if e["id"] == entry_id]

    return [_format_apa(e["type"], e["fields"]) for e in entries]


if __name__ == "__main__":
    # Example manual usage:
    #   python bibtexxml_apa.py path/to/file.xml
    import sys

    if len(sys.argv) < 2:
        raise SystemExit("Usage: python bibtexxml_apa.py <xml-file> [entry-id]")

    xml_path = sys.argv[1]
    selected_id = sys.argv[2] if len(sys.argv) > 2 else None

    for citation in bibtexxml_to_apa(xml_path, selected_id):
        print(citation)
