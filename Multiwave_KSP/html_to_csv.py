from bs4 import BeautifulSoup
import pandas as pd

input_path = "Swift_BAT_Hard_X-ray_Transient_Monitor.html"
output_path = "Swift_BAT_Transient_Sources.csv"

with open(input_path, "rb") as f:
    soup = BeautifulSoup(f.read(), "html.parser")

table = soup.find("table")
rows = table.find_all("tr")

headers = [th.get_text(strip=True) for th in rows[0].find_all(["th", "td"])]

data = []
for row in rows[1:]:
    cells = row.find_all(["td", "th"])
    row_data = [c.get_text(strip=True) for c in cells]
    if any(row_data):
        data.append(row_data)

df = pd.DataFrame(data, columns=headers)
df = df.rename(columns={"": "Row"})
df.to_csv(output_path, index=False)

print(f"Saved {len(df)} rows to {output_path}")
