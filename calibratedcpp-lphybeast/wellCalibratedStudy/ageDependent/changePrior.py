#!/usr/bin/env python3
import os
from lxml import etree as ET

# ----------------------------
# Output folder setup
# ----------------------------
input_folder = "./script"
output_folder = "./xmls"

# ----------------------------
# Helpers
# ----------------------------
def _find_age_value(cal_node, tag_name: str):
    """
    Find age value from either:
    1. Child element with tag name (e.g., <upperAge>, <lowerAge>)
    2. Parameter element with name attribute

    Returns the value attribute or text content, or None if not found.
    """
    # Try to find direct child element (e.g., <upperAge>, <lowerAge>)
    elem = cal_node.find(f"./{tag_name}")
    if elem is not None:
        value = elem.get("value")
        if value is not None:
            return value
        if elem.text:
            return elem.text.strip()

    # Fallback: look for parameter with name attribute
    for p in cal_node.findall(".//parameter"):
        if p.get("name") == tag_name:
            if p.text is None:
                return None
            return p.text.strip()
    return None


def _find_taxonset_id(cal_node):
    """
    Supports two forms:
      1. Inline child <taxa id="TaxonSetX"> — original assumption in the script
      2. Attribute reference taxa="@TaxonSetX" — what the XML actually contains

    Returns the taxonset ID string (without the leading '@'), or None if not found.
    """
    # Form 1: inline child <taxa>
    taxa = cal_node.find(".//taxa")
    if taxa is not None:
        return taxa.get("id")
    # Form 2: taxa="@TaxonSetX" attribute reference
    ref = cal_node.get("taxa")
    if ref and ref.startswith("@"):
        return ref[1:]
    return None


def _make_distribution(indent_tail: str, k: int, lower: str, upper: str, taxset: str):
    """
    Build a <distribution> element with the new MRCAPrior format:
    <distribution id="MRCAPriorX" spec="beast.base.spec.evolution.tree.MRCAPrior"
                   taxonset="@TaxonSetX" tree="@tree">
        <distr id="tX.prior" spec="beast.base.spec.inference.distribution.Uniform">
            <lower id="RealScalarParamX" spec="beast.base.spec.inference.parameter.RealScalarParam"
                   domain="Real" value="..."/>
            <upper id="RealScalarParamX" spec="beast.base.spec.inference.parameter.RealScalarParam"
                   domain="Real" value="..."/>
        </distr>
    </distribution>
    """
    dist = ET.Element(
        "distribution",
        attrib={
            "id": f"MRCAPrior{k}",
            "spec": "beast.base.spec.evolution.tree.MRCAPrior",
            "taxonset": f"@{taxset}",
            "tree": "@tree",
        },
    )

    distr = ET.SubElement(
        dist,
        "distr",
        attrib={
            "id": f"t{k}.prior",
            "spec": "beast.base.spec.inference.distribution.Uniform",
        },
    )

    lower_elem = ET.SubElement(
        distr,
        "lower",
        attrib={
            "id": f"RealScalarParam{k*10}",
            "spec": "beast.base.spec.inference.parameter.RealScalarParam",
            "domain": "Real",
            "value": str(lower),
        },
    )

    upper_elem = ET.SubElement(
        distr,
        "upper",
        attrib={
            "id": f"RealScalarParam{k*10+1}",
            "spec": "beast.base.spec.inference.parameter.RealScalarParam",
            "domain": "Real",
            "value": str(upper),
        },
    )

    dist.text = "\n" + indent_tail + "    "
    distr.text = "\n" + indent_tail + "        "
    lower_elem.tail = "\n" + indent_tail + "        "
    upper_elem.tail = "\n" + indent_tail + "    "
    distr.tail = "\n" + indent_tail

    return dist


def unwrap_calibrationprior(root):
    """
    Remove the <distribution spec="calibrationprior.CalibrationPrior"> wrapper,
    promoting its <calibration> children directly into the parent element.
    """
    count = 0
    nodes = root.xpath('.//distribution[@spec="calibrationprior.CalibrationPrior"]')
    for node in nodes:
        parent = node.getparent()
        if parent is None:
            continue
        idx = parent.index(node)
        children = list(node)
        for i, child in enumerate(children):
            parent.insert(idx + i, child)
        parent.remove(node)
        count += 1
    return count


def replace_calibrations_with_mrcapriors(root):
    """
    Replace every <calibration spec="calibrationprior.CalibrationCladePrior">
    with a <distribution spec="beast.base.spec.evolution.tree.MRCAPrior">.

    The new MRCAPrior gets attributes with the lowerAge / upperAge values,
    and a <distr> child with <lower> and <upper> elements.

    Handles both inline <taxa> children and taxa="@TaxonSetX" attribute references.
    Handles both nested <upperAge>/<lowerAge> elements and parameter-style elements.
    """
    cal_nodes = root.xpath('.//calibration[@spec="calibrationprior.CalibrationCladePrior"]')

    replaced = 0
    k = 2

    for cal in cal_nodes:
        # Look for upperAge and lowerAge as child elements (with @value) or parameters
        lower = _find_age_value(cal, "lowerAge")
        upper = _find_age_value(cal, "upperAge")
        taxset_id = _find_taxonset_id(cal)

        if lower is None or upper is None or taxset_id is None:
            print(
                f"  WARNING: skipping calibration id='{cal.get('id')}' — "
                f"missing lower={lower}, upper={upper}, taxset_id={taxset_id}"
            )
            continue

        parent = cal.getparent()
        if parent is None:
            continue

        # Determine indentation from surrounding whitespace
        indent_tail = ""
        if cal.tail and "\n" in cal.tail:
            indent_tail = cal.tail.split("\n")[-1]
        elif parent.text and "\n" in parent.text:
            indent_tail = parent.text.split("\n")[-1]

        k += 1
        new_dist = _make_distribution(indent_tail, k, lower, upper, taxset_id)
        # Preserve the original element's trailing whitespace
        new_dist.tail = cal.tail

        idx = parent.index(cal)
        parent.remove(cal)
        parent.insert(idx, new_dist)

        replaced += 1

    return replaced


def replace_filename_paths(root):
    """
    Replace fileName=" with fileName="../calibratedcpp-lphybeast/wellCalibratedStudy/ageDependent/xmls/
    in all attributes throughout the XML.
    """
    count = 0
    for elem in root.iter():
        if "fileName" in elem.attrib:
            old_value = elem.get("fileName")
            if old_value:
                # Replace fileName=" with the new path
                new_value = old_value.replace(
                    'fileName="',
                    'fileName="../calibratedcpp-lphybeast/wellCalibratedStudy/ageDependent/xmls/',
                    1
                )
                # More robust: if the attribute just contains a filename, prepend the path
                if not old_value.startswith(".."):
                    new_value = "../calibratedcpp-lphybeast/wellCalibratedStudy/ageDependent/xmls/" + old_value
                    elem.set("fileName", new_value)
                    count += 1
    return count


# ----------------------------
# Main processing
# ----------------------------
def modify_xml_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(input_folder):
        if not filename.endswith(".xml"):
            continue

        filepath = os.path.join(input_folder, filename)
        try:
            xml = ET.parse(filepath)
        except ET.XMLSyntaxError as e:
            print(f"ERROR parsing {filepath}: {e}")
            continue

        root = xml.getroot()

        # Step 1: Remove the outer CalibrationPrior distribution wrapper,
        #         promoting its <calibration> children into <prior>.
        n_unwrapped = unwrap_calibrationprior(root)

        # Step 2: Replace each <calibration CalibrationCladePrior> with
        #         a <distribution MRCAPrior> using the extracted age bounds.
        n_replaced = replace_calibrations_with_mrcapriors(root)

        # Step 3: Replace fileName paths
        n_paths = replace_filename_paths(root)

        output_path = os.path.join(output_folder, filename)
        xml.write(output_path, pretty_print=True, xml_declaration=True, encoding="UTF-8")
        print(
            f"Processed '{filename}': "
            f"{n_unwrapped} wrapper(s) unwrapped, "
            f"{n_replaced} calibration(s) replaced, "
            f"{n_paths} path(s) updated → {output_path}"
        )


# ----------------------------
# Run processing
# ----------------------------
if __name__ == "__main__":
    if not os.path.isdir(input_folder):
        print(f"Error: input folder not found: {input_folder}")
        exit(1)
    modify_xml_files(input_folder=input_folder, output_folder=output_folder)
    print(f"\nAll files processed. Output saved to: {output_folder}")