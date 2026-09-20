#!/bin/bash
#SBATCH --job-name=FABDEM


out_dir="/home/"
mkdir -p "$out_dir"
cd "$out_dir" || exit 1

python3 - <<'PY'
import re
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin

xmin, ymin, xmax, ymax = -120, -35, -30, 35 #neotropics
dataset_page = "https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn"

html = requests.get(dataset_page).text
soup = BeautifulSoup(html, "html.parser")

pattern = re.compile(r"[NS]\d{2}[EW]\d{3}-[NS]\d{2}[EW]\d{3}_FABDEM_V1-2\.zip")

def coord_to_value(coord):
    hemi = coord[0]
    value = int(coord[1:])
    return -value if hemi in ["S", "W"] else value

def tile_bounds(filename):
    tile = filename.replace("_FABDEM_V1-2.zip", "")
    a, b = tile.split("-")

    lat1 = coord_to_value(a[0:3])
    lon1 = coord_to_value(a[3:7])
    lat2 = coord_to_value(b[0:3])
    lon2 = coord_to_value(b[3:7])

    west = min(lon1, lon2)
    east = max(lon1, lon2)
    south = min(lat1, lat2)
    north = max(lat1, lat2)

    return west, south, east, north

def intersects_neotropics(bounds):
    west, south, east, north = bounds

    return not (
        east <= xmin or
        west >= xmax or
        north <= ymin or
        south >= ymax
    )

available = {}

for a in soup.find_all("a", href=True):
    text = a.get_text(" ", strip=True)
    href = a["href"]

    match = pattern.search(text) or pattern.search(href)

    if match:
        filename = match.group(0)
        url = urljoin("https://data.bris.ac.uk", href)
        available[filename] = url

selected = {
    filename: url
    for filename, url in available.items()
    if intersects_neotropics(tile_bounds(filename))
}

print("Available FABDEM files found:", len(available))
print("FABDEM tiles selected for Neotropics:", len(selected))

with open("fabdem_neotropics_urls.txt", "w") as f:
    for filename, url in sorted(selected.items()):
        f.write(url + "\n")
PY

echo "Starting parallel FABDEM download..."

cat fabdem_neotropics_urls.txt | xargs -n 1 -P 12 wget -c

echo "Testing downloaded ZIP files..."

for f in *.zip; do
    echo "Testing $f"
    unzip -t "$f" > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        echo "OK: $f"
    else
        echo "CORRUPT: $f"
    fi
done

echo "FABDEM download and ZIP check finished."
