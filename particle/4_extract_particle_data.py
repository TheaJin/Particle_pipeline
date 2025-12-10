import os
import pandas as pd
import re  # Import regex library

# Base directory containing all subjects
base_dir = "/hpc/gjin250/refine_1d/particle_development/output/"

# List to store extracted data
data = []

print("🔍 Scanning directory:", base_dir)

# Regex pattern to extract deposition values
deposition_pattern = r"TOTAL DE = ([\d.]+)\s+BRONCHIAL DE = ([\d.]+)\s+ALVEOLAR DE= ([\d.]+)\s+EXTRATHORACIC DE = ([\d.]+)"

# Loop through all directories in the base directory
for subject_folder in os.listdir(base_dir):
    if subject_folder.endswith("-cluster-12"):  # Only process folders ending with '-scaled'
        subject = subject_folder.replace("-cluster-12", "")  # Extract subject name (e.g., '009-scaled' -> '009')
        subject_path = os.path.join(base_dir, subject_folder)

        print(f"📂 Found subject: {subject_folder}")

        # Loop through all particle size folders
        for particle_size in os.listdir(subject_path):
            particle_path = os.path.join(subject_path, particle_size, "final_results.txt")

            # Check if final_results.txt exists in the particle size folder
            if os.path.exists(particle_path):
                print(f"✅ Found file: {particle_path}")

                with open(particle_path, "r") as file:
                    content = file.read()  # Read entire file content

                # Print file content for debugging
                print(f"📜 Content of {particle_path}:\n{content}")

                # Try to find the deposition values using regex
                match = re.search(deposition_pattern, content)
                if match:
                    TDE = float(match.group(1))  # Extract TOTAL DE
                    BDE = float(match.group(2))  # Extract BRONCHIAL DE
                    ADF = float(match.group(3))  # Extract ALVEOLAR DE
                    DE_ET = float(match.group(4))  # Extract EXTRATHORACIC DE

                    intra_thoracic_dep = TDE - DE_ET  # Compute intra-thoracic deposition

                    # Append to data list
                    data.append({
                        "Subject": subject,
                        "Particle Size (µm)": float(particle_size),
                        "TDE": TDE,
                        "DE_ET": DE_ET,
                        "BDE": BDE,
                        "ADF": ADF,
                        "Flow (ml)": 100,
                        "Intra-thoracic Deposition": intra_thoracic_dep
                    })

                    print(f"✅ Extracted values for {subject} at {particle_size}µm: TDE={TDE}, DE_ET={DE_ET}, BDE={BDE}, ADF={ADF}")
                else:
                    print(f"⚠️ No deposition values found in {particle_path}")

# Convert to DataFrame
df = pd.DataFrame(data)

# Save to Excel if there's data
output_file = "particle_deposition_results_cluster_12.xlsx"
if not df.empty:
    df.to_excel(output_file, index=False)
    print(f"✅ Extraction complete! Results saved to {output_file}")
else:
    print("⚠️ No data extracted. Check your folder structure and `final_results.txt` files.")
