package ahrd.view;

import static ahrd.controller.Settings.getSettings;

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
	public static final DecimalFormat FRMT = new DecimalFormat(
			"#,######0.######");

	private BufferedWriter pathBufWrtr;
	private BufferedWriter outBufWrtr;
	private List<String> sortedBlastDatabases;

	public TrainerOutputWriter() throws IOException {
		super();
		// Ensure Blast-Database-Parameters always appear in the right columns:
		this.sortedBlastDatabases = new ArrayList<String>(getSettings()
				.getBlastDatabases());
		Collections.sort(this.sortedBlastDatabases);
		// Prepare buffered output-writer:
		this.pathBufWrtr = new BufferedWriter(new FileWriter(getSettings()
				.getPathToSimulatedAnnealingPathLog()));
		// And write the header into the path-log:
		this.pathBufWrtr.write(generateHeader(false));
	}

	public String generateHeader(boolean addAvgMaxEvaluationScoreCol) {
		String hdr = "Temperature\t";
		if (addAvgMaxEvaluationScoreCol)
			hdr += "Average Maximum-Evaluation-Score\t";
		hdr += "Average Evaluation-Score\tAccepted\tAverage True-Positive-Rate\tAverage False-Positive-Rate\tDescription-Score-Pattern-Factor-Weight\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight";
		for (String blastDb : this.sortedBlastDatabases) {
			hdr += "\t" + blastDb + "-Weight";
			hdr += "\t" + blastDb + "-Description-Score-Bit-Score-Weight";
		}
		hdr += "\n";
		return hdr;
	}

	public void writeIterationOutput(Settings currentSettings, int accepted)
			throws IOException {
		this.pathBufWrtr.write(settingsRow(currentSettings, accepted, null));
	}

	/**
	 * Writes out the final output and cleanes up both used buffered Writer.
	 * 
	 * @param acceptedSettings
	 * @param avgMaxEvaluationScore
	 * @throws IOException
	 */
	public void writeFinalOutput(Settings acceptedSettings,
			Double avgMaxEvaluationScore) throws IOException {
		// Clean up buffered Sim-Anneal-Path-Log-Writer:
		this.pathBufWrtr.close();

		// Write output about found best performing Parameters:
		this.outBufWrtr = new BufferedWriter(new FileWriter(getSettings()
				.getPathToOutput()));
		// this.outBufWrtr.write("Found best scoring Parameters:\n");
		this.outBufWrtr.write(generateHeader(true));
		this.outBufWrtr.write(settingsRow(acceptedSettings, 1,
				avgMaxEvaluationScore));

		// Clean buffered Output-Writer:
		this.outBufWrtr.close();
	}

	public String settingsRow(Settings s, int accepted, Double avgMaxEvalScore) {
		String col = s.getTemperature().toString() + "\t";
		if (avgMaxEvalScore != null)
			col += avgMaxEvalScore + "\t";
		col += s.getAvgEvaluationScore() + "\t"
				+ FRMT.format(s.getAvgTruePositivesRate()) + "\t" + accepted
				+ "\t" + FRMT.format(s.getAvgFalsePositivesRate()) + "\t"
				+ FRMT.format(s.getDescriptionScorePatternFactorWeight())
				+ "\t" + FRMT.format(s.getTokenScoreBitScoreWeight()) + "\t"
				+ FRMT.format(s.getTokenScoreDatabaseScoreWeight()) + "\t"
				+ FRMT.format(s.getTokenScoreOverlapScoreWeight());
		for (String blastDb : this.sortedBlastDatabases) {
			col += "\t" + FRMT.format(s.getBlastDbWeight(blastDb));
			col += "\t"
					+ FRMT.format(s.getDescriptionScoreBitScoreWeight(blastDb));
		}
		col += "\n";
		return col;
	}
}
