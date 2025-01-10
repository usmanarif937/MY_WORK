import ROOT

def compare_and_modify_histograms(file1, file2, file3, output_file_name, modified_output_file_name):
    # X-axis labels based on histogram names
    x_axis_labels = {
        "pt": "Transverse Momentum [GeV]",
        "eta": "Pseudorapidity (#eta)",
        "phi": "Azimuthal Angle (#phi)"
    }

    # Step 1: Compare histograms and save to a ROOT file
    def compare_histograms(file1, file2, file3, output_file_name):
        # Open the ROOT files
        f1 = ROOT.TFile(file1)
        f2 = ROOT.TFile(file2)
        f3 = ROOT.TFile(file3)

        # Create a new ROOT file to store the canvases
        output_file = ROOT.TFile(output_file_name, "RECREATE")

        # List of histogram names to be compared
        hist_names = [
            "bq0_pt", "bq0_eta", "bq0_phi",
            "bq1_pt", "bq1_eta", "bq1_phi",
            "q0_pt", "q0_eta", "q0_phi",
            "q1_pt", "q1_eta", "q1_phi",
            "muon0_pt", "muon0_eta", "muon0_phi",
            "muon1_pt", "muon1_eta", "muon1_phi",
            "wboson0_pt", "wboson0_eta", "wboson0_phi",
            "wboson1_pt", "wboson1_eta", "wboson1_phi",
            "top_quark0_pt", "top_quark0_eta", "top_quark0_phi",
            "top_quark1_pt", "top_quark1_eta", "top_quark1_phi"
        ]

        # Create a canvas for each histogram comparison
        for hist_name in hist_names:
            canvas = ROOT.TCanvas(hist_name, hist_name, 800, 800)

            hist1 = f1.Get(hist_name)
            hist2 = f2.Get(hist_name)
            hist3 = f3.Get(hist_name)

            if hist1 and hist2 and hist3:
                hist1.SetLineColor(ROOT.kRed)
                hist2.SetLineColor(ROOT.kBlue)
                hist3.SetLineColor(ROOT.kGreen)

                # Disable the stats box
                hist1.SetStats(0)
                hist2.SetStats(0)
                hist3.SetStats(0)

                hist1.Draw("HIST")
                hist2.Draw("HIST SAME")
                hist3.Draw("HIST SAME")

                # Set axis titles based on histogram name
                for key, label in x_axis_labels.items():
                    if key in hist_name:
                        hist1.GetXaxis().SetTitle(label)
                        break

                # Add legend for clarity,
               
                legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85) #top right corner 
              #  legend = ROOT.TLegend(0.35, 0.2, 0.65, 0.3) # bottom center
                legend.SetBorderSize(0)
                legend.SetTextFont(43)
                legend.SetTextSize(28)
                legend.AddEntry(hist1, "NP", "l")
                legend.AddEntry(hist2, "SM", "l")
                legend.AddEntry(hist3, "SM+NP", "l")
                legend.Draw()

                # Write the canvas to the output ROOT file
                canvas.Write()
                print(f"Saved comparison for histogram '{hist_name}' in ROOT file.")
            else:
                print(f"Histogram '{hist_name}' not found in one of the files.")

        # Close files
        f1.Close()
        f2.Close()
        f3.Close()
        output_file.Close()

    # Step 2: Modify histograms and save as PDFs
    def modify_histograms(input_file_name, output_file_name):
        # Open the input ROOT file
        input_file = ROOT.TFile(input_file_name, "READ")
        
        # Create a new ROOT file to store the modified canvases
        output_file = ROOT.TFile(output_file_name, "RECREATE")

        # Get the list of keys (canvases) in the input file
        keys = input_file.GetListOfKeys()

        for key in keys:
            # Get the object name from the key
            obj_name = key.GetName()

            # Get the canvas
            canvas = input_file.Get(obj_name)
            if not isinstance(canvas, ROOT.TCanvas):
                continue

            # Retrieve all histograms drawn on the canvas
            primitives = canvas.GetListOfPrimitives()

            for primitive in primitives:
                if isinstance(primitive, ROOT.TH1):  # Check if the primitive is a histogram
                    # Modify Y-axis title
                    primitive.GetYaxis().SetTitle("Number of Events")

                    # Determine x-axis label based on the histogram name
                    for key, label in x_axis_labels.items():
                        if key in primitive.GetName():
                            primitive.GetXaxis().SetTitle(label)
                            break

                    # Adjust title fonts and sizes
                    font_code = 43
                    font_size = 28
                    primitive.GetXaxis().SetTitleFont(font_code)
                    primitive.GetXaxis().SetTitleSize(font_size)
                    primitive.GetXaxis().CenterTitle()
                    primitive.GetYaxis().SetTitleFont(font_code)
                    primitive.GetYaxis().SetTitleSize(font_size)
                    primitive.GetYaxis().CenterTitle()
                    
                    # Set the histogram title and shift it downward
                    primitive.SetTitle("Four Tops in SM & 2HDM")
                    primitive.SetTitleFont(font_code)
                  

                    # Adjust Y-axis range for better appearane
                    max_value = primitive.GetMaximum()
                    primitive.SetMaximum(max_value * 1.2)
                    primitive.SetMinimum(0)

                    # Increase line width and set bold line style
                    primitive.SetLineWidth(4)
                    primitive.SetLineStyle(1)

            # Adjust canvas size and margins
            canvas.SetCanvasSize(800, 800)
            canvas.SetTopMargin(0.06)
            canvas.SetBottomMargin(0.1)
            canvas.SetLeftMargin(0.13)
            canvas.SetRightMargin(0.05)

            # Write the modified canvas to the output ROOT file
            canvas.Write()

            # Save the canvas as a PDF
            pdf_file_name = f"{obj_name}.pdf"
            canvas.SaveAs(pdf_file_name)

        # Close the files
        input_file.Close()
        output_file.Close()

        print(f"Modified histograms saved to '{output_file_name}', and individual PDFs created.")

# Run both steps
    compare_histograms(file1, file2, file3, output_file_name)
    modify_histograms(output_file_name, modified_output_file_name)

# Example usage
compare_and_modify_histograms(
    "unweighted_events_NP.root",
    "unweighted_events_SM.root",
    "unweighted_events_SM_NP.root",
    "histogram_comparison.root",
    "modified_histogram_comparison.root"
)
   

