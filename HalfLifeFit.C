void HalfLifeFitFromHistogram() {
    // Example data points (y-values only) as a list
    std::vector<double> data = {32, 16, 8, 4, 2, 1}; // Values decreasing with half-life

    // Create a histogram based on the data
    int nBins = data.size();
    double xMin = 0;
    double xMax = nBins; // Assuming each bin corresponds to an integer time step
    TH1F *hist = new TH1F("hist", "Exponential Decay with Half-Life;Time;Value", nBins, xMin, xMax);

    // Fill the histogram with the data
    for (int i = 0; i < nBins; ++i) {
        hist->SetBinContent(i + 1, data[i]); // Bin indices start at 1 in ROOT
    }

    // Define the exponential decay function for fitting
    TF1 *expFit = new TF1("expFit", "[0] * pow(2, -x/[1])", xMin, xMax);
    expFit->SetParameters(100, 1.5); // Initial guesses: a = 100, T = 1.5

    // Fit the function to the histogram
    hist->Fit(expFit, "R");

    // Draw the histogram and the fit
    TCanvas *c1 = new TCanvas("c1", "Exponential Decay with Half-Life", 800, 600);
    gStyle->SetOptFit(1);
    hist->Draw();
    expFit->Draw("same");

    // Print the fit parameters
    double a = expFit->GetParameter(0); // Parameter a
    double T = expFit->GetParameter(1); // Half-life T
    printf("Fit parameters: a = %f, half-life = %f\n", a, T);
}
