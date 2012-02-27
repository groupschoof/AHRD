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

	private BufferedWriter bufWrtr;
	private List<String> sortedBlastDatabases;

	public TrainerOutputWriter() throws IOException {
		super();
		// Ensure Blast-Database-Parameters always appear in the right columns:
		this.sortedBlastDatabases = new ArrayList<String>(getSettings()
				.getBlastDatabases());
		Collections.sort(this.sortedBlastDatabases);
		// Prepare buffered output-writer:
		this.bufWrtr = new BufferedWriter(new FileWriter(getSettings()
				.getPathToOutput()));
		// And write the header:
		writeHeader();
	}

	public void writeHeader() throws IOException {
		this.bufWrtr
				.write("Temperature\tAverage Evaluation-Score\tAverage True-Positive-Rate\tAverage False-Positive-Rate\tDescription-Score-Pattern-Factor-Weight\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight");
		for (String blastDb : this.sortedBlastDatabases) {
			this.bufWrtr.write("\t" + blastDb + "-Weight");
			this.bufWrtr.write("\t" + blastDb
					+ "-Description-Score-Bit-Score-Weight");
		}
		this.bufWrtr.write("\n");
	}

	public void writeIterationOutput(Settings currentSettings)
			throws IOException {
		this.bufWrtr.write(settingsRow(currentSettings, false));
	}

	public void writeFinalOutput(Settings acceptedSettings) throws IOException {
		this.bufWrtr.write(settingsRow(acceptedSettings, true));
		this.bufWrtr.close();
	}

	public String settingsRow(Settings s, boolean bestSettings) {
		String col = bestSettings ? "best-settings" : s.getTemperature()
				.toString();
		col += "\t" + s.getAvgEvaluationScore() + "\t"
				+ FRMT.format(s.getAvgTruePositivesRate()) + "\t"
				+ FRMT.format(s.getAvgFalsePositivesRate()) + "\t"
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
