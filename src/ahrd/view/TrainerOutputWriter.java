package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

import ahrd.controller.Parameters;

public abstract class TrainerOutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	protected static final DecimalFormat FRMT = new DecimalFormat("#,######0.######");
	protected static final String separator = "\t";

	protected BufferedWriter pathBufWrtr;
	protected BufferedWriter outBufWrtr;

	protected TrainerOutputWriter() throws IOException {
		super();
		// Prepare buffered output-writer
		this.pathBufWrtr = new BufferedWriter(new FileWriter(getSettings().getPathToTrainingPathLog()));
	}
	
	/**
	 * Generates a header string for training path output in accordance to the type of parameters used.
	 * @param isFinalOutput
	 * @param p
	 * @return header
	 */
	public abstract String generateHeader(boolean isFinalOutput, Parameters p);
	
	/**
	 * Generates a header in accordance to the Parameters object p and writes it to the path log file
	 * @param p
	 * @throws IOException
	 */
	public void writeHeader(Parameters p) throws IOException {
		this.pathBufWrtr.write(generateHeader(false, p));
	}
	
	/**
	 * Generates a line for the training path output and writes it to file.
	 * @param iteration
	 * @param currentParameters
	 * @param diff
	 * @param trainingStat
	 * @throws IOException
	 */
	public void writeIterationOutput(Integer iteration, Parameters currentParameters, double diff, String trainingStat) throws IOException {
		String col = iteration + separator
				+ currentParameters.getAvgEvaluationScore() + separator
				+ diff + separator
				+ trainingStat + separator
				+ FRMT.format(currentParameters.getAvgPrecision()) + separator
				+ FRMT.format(currentParameters.getAvgRecall()) + separator
				+ currentParameters.formatForOutput(FRMT, separator)
				+ "\n";
		this.pathBufWrtr.write(col);
		this.pathBufWrtr.flush();
	}
	
	/**
	 * Writes out the final output and cleans up both used buffered writers.
	 * 
	 * @param iteration
	 * @param avgMaxEvaluationScore 
	 * @param acceptedParameters
	 * @throws IOException
	 */
	public void writeFinalOutput(Integer iteration, Double avrgMaxEvaluationScore, Parameters acceptedParameters) throws IOException {
		// Clean up buffered writer for the training iterations
		this.pathBufWrtr.close();
		// Initialize buffered output writer for final output
		this.outBufWrtr = new BufferedWriter(new FileWriter(getSettings().getPathToOutput()));
		// Write header
		this.outBufWrtr.write(generateHeader(true, acceptedParameters));
		// Write found best performing parameters
		String col = iteration + separator
				+ avrgMaxEvaluationScore + separator
				+ acceptedParameters.getAvgEvaluationScore() + separator
				+ FRMT.format(acceptedParameters.getAvgPrecision()) + separator
				+ FRMT.format(acceptedParameters.getAvgRecall()) + separator
				+ acceptedParameters.formatForOutput(FRMT, separator)
				+ "\n";
		this.outBufWrtr.write(col);
		// Clean up buffered Output-Writer:
		this.outBufWrtr.close();
	}
	
}
