package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import ahrd.controller.AHRD;
import ahrd.model.Blast2GoAnnot;
import ahrd.model.BlastResult;
import ahrd.model.GeneOntologyResult;
import ahrd.model.InterproResult;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;

public class OutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	public static final DecimalFormat FRMT = new DecimalFormat("#,###0.###");

	private Collection<Protein> proteins;
	
	public OutputWriter(Collection<Protein> proteins) {
		setProteins(proteins);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings()
				.getPathToOutput()));

		// Column-Names:
		bw.write("# AHRD-Version " + AHRD.VERSION + "\n");
		bw.write("\n");
		bw
				.write("Protein-Accesion\tBlast-Hit-Accession\tAHRD-Quality-Code\tHuman-Readable-Description\tInterpro-ID (Description)\tGene-Ontology-ID (Name)");

		if (getSettings().isInTrainingMode()) {
			bw
					.write("\tHRD-Length\tReference-Description\tRef-Lenght\tEvaluation-Score\tDiff-to-bestCompetitor\tTPR\tFPR");
		}
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			bw.write(buildBestBlastHitsHeader());
		}
		if (getSettings().getWriteTokenSetToOutput()) {
			bw.write("\t\"Tokens (tkn->score)\"");
		}
		if (getSettings().getWriteScoresToOutput()) {
			bw
					.write("\tSum(Token-Scores)\tTokenHighScore\tCorrection-Factor\tGO-Score\tLexical-Score\tRelativeBitScore\tDescriptionLineFrequency\tMax(DescLineFreq)\tPattern-Factor");
		}
		if (getSettings().getPathToBlast2GoAnnotations() != null
				&& !getSettings().getPathToBlast2GoAnnotations().equals("")) {
			bw
					.write("\tBlast2GO-Annotation\tBlast2GO-Length\tBlast2GO-Evaluation-Score");
		}

		bw.write("\n");

		for (Protein prot : getProteins()) {
			// Generate the Human Readable Description:
			String csvRow = buildDescriptionLine(prot);

			// If in Evaluator-Mode write out the Evaluator-Score and the
			// Reference-Description:
			if (getSettings().isInTrainingMode()) {
				csvRow += buildTrainerColumns(prot);
			}
			// Append further information, if requested:
			if (getSettings().getWriteBestBlastHitsToOutput()) {
				csvRow += buildBestBlastHitsColumns(prot);
			}
			if (getSettings().getWriteTokenSetToOutput()) {
				csvRow += buildTokenSetCell(prot);
			}
			if (getSettings().getWriteScoresToOutput()) {
				csvRow += buildDescScoreCells(prot);
			}
			if (getSettings().getPathToBlast2GoAnnotations() != null
					&& !getSettings().getPathToBlast2GoAnnotations().equals("")) {
				csvRow += buildBlast2GoColumns(prot);
			}

			// Write row to CSV:
			csvRow += "\n";
			bw.write(csvRow);
		}

		bw.close();
	}

	public String buildBlast2GoColumns(Protein prot) {
		String csvCols = "";
		List<Blast2GoAnnot> rankedBlast2GoAnnots = prot
				.getEvaluationScoreCalculator().sortBlast2GoAnnotsByEvalScore();
		if (rankedBlast2GoAnnots != null && !rankedBlast2GoAnnots.isEmpty()) {
			Blast2GoAnnot bestB2ga = rankedBlast2GoAnnots
					.get(rankedBlast2GoAnnots.size() - 1);
			csvCols += "\t" + bestB2ga.getDescription() + "\t"
					+ bestB2ga.getEvaluationTokens().size() + "\t"
					+ FRMT.format(bestB2ga.getEvaluationScore());
		} else {
			csvCols += "\t\t0\t0.0";
		}
		return csvCols;
	}

	/**
	 * @param Protein
	 *            prot
	 * @return String - Part of the CSV-Row with the columns Evaluator-Score and
	 *         Reference-Description.
	 */
	public String buildTrainerColumns(Protein prot) {
		// HEADER:
		// \tHRD-Length\tReference-Description\tRef-Lenght\tEvaluation-Score\tDiff-to-bestCompetitor
		String csvCells = "";
		// HRD-Length reference and AHRD's performance:
		if (prot.getEvaluationScoreCalculator().getEvalutionScore() != null) {
			// HRD-Length ref f1score diff-to-best-competitor:
			csvCells += "\t";
			if (prot.getDescriptionScoreCalculator()
					.getHighestScoringBlastResult() != null)
				csvCells += prot.getDescriptionScoreCalculator()
						.getHighestScoringBlastResult().getEvaluationTokens()
						.size();
			else
				csvCells += "0";
			csvCells += "\t"
					+ prot.getEvaluationScoreCalculator()
							.getReferenceDescription().getDescription()
					+ "\t"
					+ prot.getEvaluationScoreCalculator()
							.getReferenceDescription().getTokens().size()
					+ "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator()
							.getEvalutionScore())
					+ "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator()
							.getEvalScoreMinBestCompScore())
					+ "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator()
							.getTruePositivesRate())
					+ "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator()
							.getFalsePositivesRate());
		} else
			csvCells = "\t\t\t\t\t\t\t";
		return csvCells;
	}

	/**
	 * Append the following columns to the current CSV-Row:
	 * sum_of_all_token_scores, token_high_score, correction_factor, go_score,
	 * lexical_score, rel_bit_score, desc_line_frequency,
	 * max_desc_line_frequency, pattern_factor
	 * 
	 * @note: Printing out these values may slow down the running-time, as these
	 *        values are not necessarily stored in memory.
	 * 
	 * @param prot
	 * @return String - Part of the CSV-Row holding the values of the above
	 *         columns.
	 */
	public String buildDescScoreCells(Protein prot) {
		String csvCells = "";
		// Found a high scoring description?
		if (prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() == null) {
			csvCells = "\t\t\t\t\t\t\t\t";
		} else {
			BlastResult hsbr = prot.getDescriptionScoreCalculator()
					.getHighestScoringBlastResult();
			csvCells += "\t"
					+ FRMT.format(prot.getTokenScoreCalculator()
							.sumOfAllTokenScores(hsbr));
			csvCells += "\t"
					+ FRMT.format(prot.getTokenScoreCalculator()
							.getTokenHighScore());
			csvCells += "\t"
					+ FRMT.format(prot.getLexicalScoreCalculator()
							.correctionFactor(hsbr));
			csvCells += "\t"
					+ FRMT.format(prot.getLexicalScoreCalculator()
							.geneOntologyScore(hsbr));
			csvCells += "\t"
					+ FRMT.format(prot.getLexicalScoreCalculator()
							.lexicalScore(hsbr));
			csvCells += "\t"
					+ FRMT.format(prot.getDescriptionScoreCalculator()
							.relativeBlastScore(hsbr));
			csvCells += "\t"
					+ FRMT.format(prot.getDescriptionScoreCalculator()
							.getDescLinePatternFrequencies().get(
									hsbr.patternize()));
			csvCells += "\t"
					+ FRMT.format(prot.getDescriptionScoreCalculator()
							.getMaxDescriptionLineFrequency());
			csvCells += "\t"
					+ FRMT.format(prot.getDescriptionScoreCalculator()
							.patternFactor(hsbr));
		}
		return csvCells;
	}

	public String buildTokenSetCell(Protein prot) {
		String tokenSetCell = "\t";

		for (String token : prot.getTokenScoreCalculator().getTokenScores()
				.keySet()) {
			tokenSetCell += "["
					+ token
					+ "->"
					+ FRMT.format(prot.getTokenScoreCalculator()
							.getTokenScores().get(token)) + "]";
		}

		return tokenSetCell;
	}

	public String buildBestBlastHitsHeader() {
		String hdr = "";
		for (String blastDb : getSettings().getBlastDatabases()) {
			if (blastDb != null && !blastDb.equals(""))
				hdr += ("\tBest BlastHit against '" + blastDb + "'");
			if (getSettings().isInTrainingMode())
				hdr += "\tLength\tEvaluation-Score";
		}
		return hdr;
	}

	public String buildBestBlastHitsColumns(Protein prot) {
		String csvRow = "";

		for (String blastDb : getSettings().getBlastDatabases()) {
			if (prot.getEvaluationScoreCalculator().getUnchangedBlastResults()
					.get(blastDb) != null) {
				BlastResult bestBr = prot.getEvaluationScoreCalculator()
						.getUnchangedBlastResults().get(blastDb);
				csvRow += "\t\"" + bestBr.getAccession() + " "
						+ bestBr.getDescription() + "\"";
				if (bestBr.getEvaluationScore() != null)
					csvRow += "\t" + bestBr.getEvaluationTokens().size() + "\t"
							+ FRMT.format(bestBr.getEvaluationScore());
			} else {
				csvRow += "\t";
				if (getSettings().isInTrainingMode())
					csvRow += "\t0\t0.0";
			}
		}
		return csvRow;
	}

	public String buildDescriptionLine(Protein protein) {
		String descLine = protein.getAccession() + "\t";
		// Blast-Results
		if (protein.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult() != null) {
			BlastResult br = protein.getDescriptionScoreCalculator()
					.getHighestScoringBlastResult();
			descLine += br.getAccession() + "\t" + qualityCode(protein) + "\t"
					+ br.getDescription() + "\t";
		} else {
			descLine += "\t\tUnknown protein\t";
		}
		// Interpro
		List<InterproResult> sortedIprs = new ArrayList<InterproResult>(protein
				.getInterproResults());
		Collections.sort(sortedIprs);
		for (Iterator<InterproResult> i = sortedIprs.iterator(); i.hasNext();) {
			InterproResult ipr = i.next();
			descLine += ipr.getId() + " (" + ipr.getName() + ")";
			if (i.hasNext())
				descLine += ", ";
		}
		descLine += "\t";
		// Gene-Ontology-Results:
		List<GeneOntologyResult> sortedGOs = new ArrayList<GeneOntologyResult>(
				protein.getGoResults());
		Collections.sort(sortedGOs);
		for (Iterator<GeneOntologyResult> i = sortedGOs.iterator(); i.hasNext();) {
			GeneOntologyResult gor = i.next();
			descLine += gor.getAcc() + " (" + gor.getName() + ")";
			if (i.hasNext())
				descLine += ", ";
		}
		return descLine;
	}

	/**
	 * Four Positions. Each gets a '*', if...
	 * 
	 * Position One: BitScore > 50 and EValue < 0.1
	 * 
	 * Position Two: Overlap > 60%
	 * 
	 * Position Three: DescriptionScore >= 0.5
	 * 
	 * Position Four: DescriptionLine shares tokens with predicted
	 * Gene-Ontology-Terms
	 * 
	 * @return the quality code
	 */
	public String qualityCode(Protein p) {
		BlastResult hsbr = p.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult();
		String qc = "";
		// Position 1
		qc += (hsbr.getBitScore() > 50.0 && hsbr.getEValue() < 0.1) ? "*" : "-";
		// Position 2
		qc += (TokenScoreCalculator.overlapScore(hsbr.getStart(),
				hsbr.getEnd(), p.getSequenceLength()) > 0.6) ? "*" : "-";
		// Position 3
		qc += (p.getDescriptionScoreCalculator().getDescriptionHighScore() >= 0.5) ? "*"
				: "-";
		// Position 4
		qc += (p.getLexicalScoreCalculator().geneOntologyScore(hsbr) > 0.0) ? "*"
				: "-";
		// Internal DescriptionScore:
		qc += "[" + FRMT.format(hsbr.getDescriptionScore()) + "]";

		return qc;
	}

	public Collection<Protein> getProteins() {
		return proteins;
	}

	public void setProteins(Collection<Protein> proteins) {
		this.proteins = proteins;
	}
}