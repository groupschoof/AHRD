package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ahrd.controller.Settings;
import ahrd.model.Fscore;

public class PrecisionRecallCurveOutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	private final DecimalFormat FRMT = new DecimalFormat("#,######0.######");

	private BufferedWriter outBufWrtr;
	private List<String> goCentricTerms = new LinkedList<String>();

	public PrecisionRecallCurveOutputWriter() throws IOException {
		super();
		// Prepare buffered output-writer:
		this.outBufWrtr = new BufferedWriter(new FileWriter(getSettings().getPathToGoAnnotationPrecisionRecallCurveFile()));
		// Load list of GO centric terms
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoTermCentricTermsFile()) {
			for (String goCentricTermFileEntry : getSettings().getGoTermCentricTerms()) {
				Pattern p = Settings.GO_TERM_CENTRIC_TERMS_FILE_GOTERM_REGEX;
				Matcher m = p.matcher(goCentricTermFileEntry);
				if (m.find()) {
					String termAcc = m.group("goTerm");
					this.goCentricTerms.add(termAcc);
				}
			}
		}
		// And write the header into the path-log:
		this.outBufWrtr.write(generateHeader());
	}

	/**
	 * Creates the first line of the file containing the column names
	 * @return header
	 * @throws IOException 
	 */
	public String generateHeader() throws IOException {
		String header = "Threshold";
		if (getSettings().doCalculateSimpleGoF1Scores()) {
			header += "\tGO-Annotations-Simple-F-Score-Precision\tGO-Annotations-Simple-F-Score-Recall";
		}
		if (getSettings().doCalculateAncestryGoF1Scores()) {
			header += "\tGO-Annotations-Ancestry-F-Score-Precision\tGO-Annotations-Ancestry-F-Score-Recall";
		}
		if (getSettings().doCalculateSemSimGoF1Scores()) {
			header += "\tGO-Annotations-SemSim-F-Score-Precision\tGO-Annotations-SemSim-F-Score-Recall";
		}
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoTermCentricTermsFile()) {
			for (String termAcc : goCentricTerms) {
				header += "\t" + termAcc + "-Precision" + "\t" + termAcc + "-Recall";
			}
		}
		header += "\n";
		return header;
	}
	
	/**
	 * Creates and writes a line to the file. 
	 * Each line starts with the threshold and contains precision and recall values depending on the scores specified in the setting 
	 * @param threshold
	 * @param simplePrecision
	 * @param simpleRecall
	 * @param ancestryPrecision
	 * @param ancestryRecall
	 * @param semSimPrecision
	 * @param semSimRecall
	 * @throws IOException
	 */
	public void writeIterationOutput(double threshold, 
			double simplePrecision, double simpleRecall, 
			double ancestryPrecision, double ancestryRecall, 
			double semSimPrecision, double semSimRecall) throws IOException {
		String line = generateBasicIterationLine(threshold, 
				simplePrecision, simpleRecall, 
				ancestryPrecision, ancestryRecall, 
				semSimPrecision, semSimRecall);
		line += "\n";
		this.outBufWrtr.write(line); 
	}
	
	public void writeIterationOutput(double threshold, 
			double simplePrecision, double simpleRecall, 
			double ancestryPrecision, double ancestryRecall, 
			double semSimPrecision, double semSimRecall,
			Map<String, Fscore> goCentricScores) throws IOException {
		String line = generateBasicIterationLine(threshold, 
				simplePrecision, simpleRecall, 
				ancestryPrecision, ancestryRecall, 
				semSimPrecision, semSimRecall);
		for (String term : this.goCentricTerms) {
			line += "\t" + FRMT.format(goCentricScores.get(term).getPrecision());
			line += "\t" + FRMT.format(goCentricScores.get(term).getRecall());
		}
		line += "\n";
		this.outBufWrtr.write(line); 
	}
	
	private String generateBasicIterationLine(double threshold, 
			double simplePrecision, double simpleRecall, 
			double ancestryPrecision, double ancestryRecall, 
			double semSimPrecision, double semSimRecall) {
		String line = FRMT.format(threshold);
		if (getSettings().doCalculateSimpleGoF1Scores()) {
			line += "\t" + FRMT.format(simplePrecision) + "\t" + FRMT.format(simpleRecall);
		}
		if (getSettings().doCalculateAncestryGoF1Scores()) {
			line += "\t" + FRMT.format(ancestryPrecision) + "\t" + FRMT.format(ancestryRecall);
		}
		if (getSettings().doCalculateSemSimGoF1Scores()) {
			line += "\t" + FRMT.format(semSimPrecision) + "\t" + FRMT.format(semSimRecall);
		}
		return line;
	}

	/**
	 * Cleans up the buffered Writer.
	 * 
	 * @throws IOException
	 */
	public void cleanUp() throws IOException {
		this.outBufWrtr.close();
	}

}