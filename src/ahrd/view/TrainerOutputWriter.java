package ahrd.view;

import static ahrd.controller.Settings.getSettings;
import static ahrd.view.AbstractOutputWriter.formattedNumberToString;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ahrd.controller.Settings;

public class TrainerOutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	public static final DecimalFormat FRMT = new DecimalFormat("#,######0.######");

	private BufferedWriter pathBufWrtr;
	private BufferedWriter outBufWrtr;
	private List<String> sortedBlastDatabases;

	public TrainerOutputWriter() throws IOException {
		super();
		// Ensure Blast-Database-Parameters always appear in the right columns:
		this.sortedBlastDatabases = new ArrayList<String>(getSettings().getBlastDatabases());
		Collections.sort(this.sortedBlastDatabases);
		// Prepare buffered output-writer:
		this.pathBufWrtr = new BufferedWriter(new FileWriter(getSettings().getPathToSimulatedAnnealingPathLog()));
		// And write the header into the path-log:
		this.pathBufWrtr.write(generateHeader(false));
	}

	public String generateHeader(boolean isFinalOutput) {
		String hdr = "Temperature\t";
		if (isFinalOutput)
			hdr += "Average Maximum-Evaluation-Score\t";
		hdr += "Average Evaluation-Score(F-Score)";
		if (!isFinalOutput)
			hdr += "\tDiff-to-curr-Accepted\tAccepted";
		hdr += "\tAverage True-Positive-Rate\tAverage False-Positive-Rate\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight";
		for (String blastDb : this.sortedBlastDatabases) {
			hdr += "\t" + blastDb + "-Weight";
			hdr += "\t" + blastDb + "-Description-Score-Bit-Score-Weight";
		}
		hdr += "\n";
		return hdr;
	}

	public void writeIterationOutput(Settings currentSettings, double diffAvgEvalScoreToCurrAccepted, int accepted)
			throws IOException {
		this.pathBufWrtr.write(settingsRow(currentSettings, diffAvgEvalScoreToCurrAccepted, accepted));
	}

	/**
	 * Writes out the final output and cleanes up both used buffered Writer.
	 * 
	 * @param acceptedSettings
	 * @param avgMaxEvaluationScore
	 * @throws IOException
	 */
	public void writeFinalOutput(Settings acceptedSettings, Double avgMaxEvaluationScore,
			Integer acceptedSettingsFoundAtTemperature) throws IOException {
		// Clean up buffered Sim-Anneal-Path-Log-Writer:
		this.pathBufWrtr.close();

		// Write output about found best performing Parameters:
		this.outBufWrtr = new BufferedWriter(new FileWriter(getSettings().getPathToOutput()));
		// this.outBufWrtr.write("Found best scoring Parameters:\n");
		this.outBufWrtr.write(generateHeader(true));
		this.outBufWrtr
				.write(finalSettingsRow(acceptedSettings, acceptedSettingsFoundAtTemperature, avgMaxEvaluationScore));

		// Clean buffered Output-Writer:
		this.outBufWrtr.close();
	}

	public String settingsRow(Settings s, double diffAvgEvalScoreToCurrAccepted, int accepted) {
		String col = s.getTemperature().toString() + "\t" + s.getAvgEvaluationScore() + "\t"
				+ diffAvgEvalScoreToCurrAccepted + "\t" + accepted + "\t"
				+ formattedNumberToString(s.getAvgTruePositivesRate()) + "\t"
				+ formattedNumberToString(s.getAvgFalsePositivesRate()) + "\t" + "\t"
				+ formattedNumberToString(s.getTokenScoreBitScoreWeight()) + "\t"
				+ formattedNumberToString(s.getTokenScoreDatabaseScoreWeight()) + "\t"
				+ formattedNumberToString(s.getTokenScoreOverlapScoreWeight());
		for (String blastDb : this.sortedBlastDatabases) {
			col += "\t" + FRMT.format(s.getBlastDbWeight(blastDb));
			col += "\t" + FRMT.format(s.getDescriptionScoreBitScoreWeight(blastDb));
		}
		col += "\n";
		return col;
	}

	public String finalSettingsRow(Settings s, Integer sFoundAtTemp, Double avgMaxEvalScore) {
		String col = sFoundAtTemp + "\t" + avgMaxEvalScore + "\t" + s.getAvgEvaluationScore() + "\t"
				+ formattedNumberToString(s.getAvgTruePositivesRate()) + "\t"
				+ formattedNumberToString(s.getAvgFalsePositivesRate()) + "\t"
				+ formattedNumberToString(s.getTokenScoreBitScoreWeight()) + "\t"
				+ formattedNumberToString(s.getTokenScoreDatabaseScoreWeight()) + "\t"
				+ formattedNumberToString(s.getTokenScoreOverlapScoreWeight());
		for (String blastDb : this.sortedBlastDatabases) {
			col += "\t" + FRMT.format(s.getBlastDbWeight(blastDb));
			col += "\t" + FRMT.format(s.getDescriptionScoreBitScoreWeight(blastDb));
		}
		col += "\n";
		return col;
	}
}
