import pandas as pd

# Load the original TSV
df = pd.read_csv("busco_results.tsv", sep="\t")

# Add new columns with default placeholder values (e.g., NaN or 0)
df["unbinned_percentage"] = None
df["absolute_unbinned_count"] = None

# Save the updated TSV
df.to_csv("unbinned_res.tsv", sep="\t", index=False)
