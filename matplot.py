import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os

# Define ROOT files and corresponding labels
root_files = {
    "unweighted_events_NP.root": "NP",
    "unweighted_events_SM.root": "SM",
    "unweighted_events_SM_NP.root": "SM+NP"
}

# Define custom colors and line styles
colors = {
    "NP": "#FFA500",  # Orange
    "SM": "#87CEEB",  # Sky Blue
    "SM+NP": "#90EE90"  # Light Green
}

line_styles = {
    "NP": "-",  # Solid line
    "SM": "--",  # Long dashed line
    "SM+NP": ":"  # Short dashed line
}

# Create output directory if it doesn't exist
output_dir = "output_plots"
os.makedirs(output_dir, exist_ok=True)

# Function to list all histograms in a ROOT file
def list_histograms(root_file):
    file = ROOT.TFile(root_file, "READ")
    hist_names = [key.GetName() for key in file.GetListOfKeys() if isinstance(file.Get(key.GetName()), ROOT.TH1)]
    file.Close()
    return hist_names

# Get the list of histograms from all files
all_hist_names = set()
for file in root_files.keys():
    all_hist_names.update(list_histograms(file))

# Function to extract histogram data from a ROOT file
def get_hist_data(root_file, hist_name):
    file = ROOT.TFile(root_file, "READ")
    hist = file.Get(hist_name)

    if not hist or not isinstance(hist, ROOT.TH1):
        file.Close()
        return None, None

    bins = np.array([hist.GetBinCenter(i) for i in range(1, hist.GetNbinsX() + 1)])
    values = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])

    file.Close()
    return bins, values

# Function to set x-axis title dynamically
def get_x_label(hist_name):
    hist_name_lower = hist_name.lower()
    if "phi" in hist_name_lower:
        return r"$\phi$"
    elif "eta" in hist_name_lower:
        return r"$\eta$"
    elif "pt" in hist_name_lower:
        return r"$P_{T}$"
    else:
        return "Bin Center"

# Loop through all histograms and plot each one
for hist_name in all_hist_names:
    plt.figure(figsize=(8, 8))  # Square canvas

    found_data = False  # Flag to check if at least one valid histogram is found

    for root_file, label in root_files.items():
        bins, values = get_hist_data(root_file, hist_name)
        if bins is not None and values is not None:
            found_data = True
            plt.step(
                bins, values, label=f"{label}",
                linestyle=line_styles.get(label, "-"),  # Apply custom line styles
                where="mid", linewidth=3,
                color=colors.get(label, "black")  # Apply custom colors
            )

    if found_data:
        # Set axis labels dynamically
        plt.xlabel(get_x_label(hist_name), fontsize=16)
        plt.ylabel("Number of Events", fontsize=16)

        # Set the main canvas title
        plt.title("Four Tops in SM and 2HDM", fontsize=18, fontweight="bold")

        # Increase legend size and box
        legend = plt.legend(loc="best", fontsize=15, frameon=True)
        legend.get_frame().set_linewidth(1.5)  # Thicker legend box

        # Set full black border around the plot
        plt.gca().spines["top"].set_color("black")
        plt.gca().spines["bottom"].set_color("black")
        plt.gca().spines["left"].set_color("black")
        plt.gca().spines["right"].set_color("black")

        # Set font size for x-axis and y-axis numbers
        plt.xticks(fontsize=14)  # Adjust x-axis number size
        plt.yticks(fontsize=14)  # Adjust y-axis number size

        # Save plot
        output_path = os.path.join(output_dir, f"{hist_name}_comparison.png")
        plt.savefig(output_path)
        print(f"✅ Saved plot: {output_path}")

        plt.close()  # Close the figure to free memory
    else:
        print(f"❌ Skipping '{hist_name}' (not found in any files)")

print("✅ All histograms have been plotted and saved!")

