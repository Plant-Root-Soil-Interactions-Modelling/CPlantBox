import csv
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.resolve().parent

def clean_value(value: str) -> str:
    """Replace '999,00' with empty string and convert commas to dots."""
    if not value:
        return ""
    v = value.strip()
    if v == "999,00":
        return ""
    return v.replace(",", ".")

def normalize_header(header: str) -> str:
    """Normalize header by stripping whitespace and removing \r and \n."""
    return header.strip().replace("\r", "").replace("\n", "").replace(" ", "").lower()

def get_out_headers(fieldnames: list[str]) -> tuple[dict[str, str], list[str]]:
    """
    Build a mapping from CSV input headers to output headers
    and return a list of output headers in order.
    """
    key_map = {}
    soil_keys = []

    for h in fieldnames:
        nh = normalize_header(h)

        if nh == "tag" or nh == "tagab1.1.":
            key_map[h] = "day_of_year"
        elif nh == "tagimmonat":
            key_map[h] = "day"
        elif nh == "monat":
            key_map[h] = "month"
        elif nh == "woche":
            key_map[h] = "woche"
        elif nh.startswith("sndschlag"):
            key_map[h] = "precipitation"
        elif nh.startswith("snettostra"):
            key_map[h] = "radiation_net"
        elif nh.startswith("sglobalstra"):
            key_map[h] = "radiation_global"
        elif nh.startswith("lufttemp["):
            key_map[h] = "air_temp_average"
        elif nh.startswith("luftfeuchte"):
            key_map[h] = "humidity"
        elif nh.startswith("windgesch"):
            key_map[h] = "wind_speed"
        elif nh.startswith("bodentemp"):
            soil_name = f"soil_temp_{nh[-4:]}"
            key_map[h] = soil_name
            soil_keys.append(h)
        elif nh.startswith("lufttempmax"):
            key_map[h] = "air_temp_max"
        elif nh.startswith("lufttempmin"):
            key_map[h] = "air_temp_min"
        elif nh.startswith("sndschlagab7uhr"):
            key_map[h] = "precipitation_since_7am"

    out_fieldnames = []
    for h in fieldnames:
        if h in key_map:
            val = key_map[h]
            if val not in out_fieldnames:
                out_fieldnames.append(val)

    return key_map, out_fieldnames

def row_to_output(row: dict[str, str], key_map: dict[str, str]) -> dict[str, str]:
    """
    Convert a single input CSV row dict to the cleaned output dict.
    """
    out_row = {}
    for in_key, out_key in key_map.items():
        value = row.get(in_key, "")
        out_row[out_key] = clean_value(value)
    return out_row

def process_csv(in_path: Path, out_path: Path) -> None:
    """Process the input CSV and write the cleaned output CSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        in_path.open("r", encoding="cp1252", newline="") as f,
        out_path.open("w", encoding="utf-8", newline="") as out_f,
    ):
        reader = csv.DictReader(f, delimiter=";")
        key_map, out_fieldnames = get_out_headers(reader.fieldnames)

        writer = csv.DictWriter(out_f, fieldnames=out_fieldnames)
        writer.writeheader()

        for row in reader:
            out_row = row_to_output(row, key_map)

            writer.writerow(out_row)

def create_clean_climate_csv(year: int) -> None:
    csv_in_path = Path(f"{BASE_DIR}/data_processed/Ruthe_{year}.CSV")
    csv_out_path = Path(
        f"{BASE_DIR}/data_processed/climate_data_{year}.csv"
    )
    csv_out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Processing soil temperatures for year {year}...")
    print(f"Input: {csv_in_path}")
    print(f"Output: {csv_out_path}")

    process_csv(csv_in_path, csv_out_path)

def main() -> None:
    for year in range(94, 99):
        create_clean_climate_csv(year)


if __name__ == "__main__":
    main()
            
            