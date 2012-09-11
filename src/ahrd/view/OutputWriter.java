package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import ahrd.controller.AHRD;
import ahrd.model.Blast2GoAnnot;
import ahrd.model.BlastResult;
import ahrd.model.Protein;

public class OutputWriter extends AbstractOutputWriter {

	protected BufferedWriter hrdScoresWriter;

	public OutputWriter(Collection<Protein> proteins) {
		super(proteins);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings()
				.getPathToOutput()));
		if (getSettings().doWriteHRDScoresToOutput())
			writeHRDScoresOutputHeader();

		// Column-Names:
		bw.write("# AHRD-Version " + AHRD.VERSION + "\n");
		bw.write("\n");
		bw
				.write("Protein-Accession-1\tBlast-Hit-Accession-2\tAHRD-Quality-Code-3\tHuman-Readable-Description-4\tInterpro-ID (Description)-5\tGene-Ontology-ID (Name)-6");

		if (getSettings().isInTrainingMode()) {
			bw
					.write("\tHRD-Length-7\tReference-Description-8\tRef-Length-9\tEvaluation-Score-10\tDiff-to-bestCompetitor-11\tTPR-12\tFPR-13");
		}
		if (getSettings().isWriteDomainArchitectureSimilarityScoresToOutput()) {
			bw.write("\tProtein-Domain-Weight-Vector-14");
		}
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			bw.write(buildBestBlastHitsHeader());
		}
		if (getSettings().getWriteTokenSetToOutput()) {
			bw.write("\t\"Tokens (tkn->score)-30\"");
		}
		if (getSettings().getWriteScoresToOutput()) {
			bw
					.write("\tSum(Token-Scores)-31\tTokenHighScore-32\tCorrection-Factor-33\tGO-Score-34\tLexical-Score-35\tRelativeBitScore-36\tDescriptionLineFrequency-37\tMax(DescLineFreq)-38\tPattern-Factor-39");
		}
		if (getSettings().getPathToBlast2GoAnnotations() != null
				&& !getSettings().getPathToBlast2GoAnnotations().equals("")) {
			bw
					.write("\tBlast2GO-Annotation\tBlast2GO-Length\tBlast2GO-Evaluation-Score");
		}
		if (getSettings().doFindHighestPossibleEvaluationScore()) {
			bw.write("\tHighest-Blast-Hit-Evaluation-Score");
		}

		bw.write("\n");

		for (Protein prot : getProteins()) {
			// Generate the Human Readable Description:
			String csvRow = buildDescriptionLine(prot, "\t");

			// If in Evaluator-Mode write out the Evaluator-Score and the
			// Reference-Description:
			if (getSettings().isInTrainingMode()) {
				csvRow += buildTrainerColumns(prot);
			}
			// Append further information, if requested:

			if (getSettings()
					.isWriteDomainArchitectureSimilarityScoresToOutput()) {
				csvRow += buildDomainWeightColumn(prot);
			}

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
			if (getSettings().doFindHighestPossibleEvaluationScore()) {
				csvRow += buildHighestPossibleEvaluationScoreColumn(prot);
			}

			// Write row to CSV:
			csvRow += "\n";
			bw.write(csvRow);

			// If AHRD is requested to write out the AHRD-Score of each
			// BlastHit's Description, do so into another file:
			if (getSettings().doWriteHRDScoresToOutput())
				writeHrdScoresOutput(prot);
		}

		// CLEAN UP:
		bw.close();
		if (getSettings().doWriteHRDScoresToOutput())
			this.hrdScoresWriter.close();
	}

	/**
	 * AHRD can be requested to log all final AHRD-Scores of each BlastHit's
	 * Description. This is needed for fitting a generalized extreme value
	 * distribution to the AHRD-Scores. The fitted gevd can later be used to
	 * calculate P-Values for the AHRD-Scores assigned to each BlastHit's
	 * Description. This method initializes the OutputWriter for the above
	 * scores.
	 */
	public void writeHRDScoresOutputHeader() throws IOException {
		// Initialize OutputWriter:
		hrdScoresWriter = new BufferedWriter(new FileWriter(getSettings()
				.getPathToHRDScoresOutput()));
		hrdScoresWriter
				.write("Protein-Accession\tBlast-Hit-Accession\tAHRD-Score\n");
	}

	/**
	 * AHRD can be requested to log all final AHRD-Scores of each BlastHit's
	 * Description. This is needed for fitting a generalized extreme value
	 * distribution to the AHRD-Scores. The fitted gevd can later be used to
	 * calculate P-Values for the AHRD-Scores assigned to each BlastHit's
	 * Description. This method writes the above scores for the argument
	 * protein.
	 * 
	 * @throws IOException
	 */
	public void writeHrdScoresOutput(Protein prot) throws IOException {
		for (String blastDatabaseName : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDatabaseName)) {
				this.hrdScoresWriter.write(prot.getAccession() + "\t"
						+ br.getAccession() + "\t" + br.getDescriptionScore()
						+ "\n");
			}
		}
	}

	public String buildHighestPossibleEvaluationScoreColumn(Protein prot) {
		return "\t"
				+ FRMT.format(prot.getEvaluationScoreCalculator()
						.getHighestPossibleEvaluationScore());
	}

	public String buildDomainWeightColumn(Protein prot) {
		String dwc = "\t";
		if (prot.getDomainWeights() != null
				&& !prot.getDomainWeights().isEmpty())
			dwc += prot.getDomainWeights().toString();
		else
			dwc += "NA-14";
		return dwc;
	}

	public String buildDomainWeightColumn(BlastResult br) {
		String dwc = "\t";
		if (br != null && br.getDomainWeights() != null
				&& !br.getDomainWeights().isEmpty())
			dwc += br.getDomainWeights().toString();
		else
			dwc += "NA-bdw";
		dwc += "\t";
		if (br != null && br.getDomainSimilarityScore() != null)
			dwc += FRMT.format(br.getDomainSimilarityScore());
		else
			dwc += "NA-bdss";
		//dwc += "\t";
		return dwc;
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
		// \tHRD-Length\tReference-Description\tRef-Length\tEvaluation-Score\tDiff-to-bestCompetitor
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
			csvCells = "\t\t8\t9\t10\t11\t12\t13";
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
			csvCells = "\t32\t33\t34\t35\t36\t37\t38\t39";
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
				hdr += ("\tBest BlastHit (BH) against " + blastDb);
			if (getSettings().isInTrainingMode())
				hdr += ("\tLength-" + blastDb)
						+ ("\tEvaluation-Score-" + blastDb);
			if (getSettings().isToComputeDomainSimilarities()) {
				hdr += ("\tBH-Domain-Weight-Vector-" + blastDb)
						+ ("\tBH-Domain-Architecture-Similarity-Score-" + blastDb);
			}
		}
		return hdr;
	}

	public String buildBestBlastHitsColumns(Protein prot) {
		String csvRow = "";
		for (String blastDb : getSettings().getBlastDatabases()) {
			BlastResult bestBr = null;
			if (prot.getEvaluationScoreCalculator().getUnchangedBlastResults()
					.get(blastDb) != null) {
				bestBr = prot.getEvaluationScoreCalculator()
						.getUnchangedBlastResults().get(blastDb);
				csvRow += "\t\"" + bestBr.getAccession() + " "
						+ bestBr.getDescription() + "\"";
				if (bestBr.getEvaluationScore() != null) {
					csvRow += "\t" + bestBr.getEvaluationTokens().size() + "\t"
							+ FRMT.format(bestBr.getEvaluationScore());
				} else
					System.out.println("MISSING Eval-Score for '"
							+ bestBr.getAccession() + "'");
			} else {
				csvRow += "\tblhit";
				if (getSettings().isInEvaluationMode()) {
					csvRow += "\t";
				}
				if (getSettings().isInTrainingMode()) {
					csvRow += "\tbllth0\tblevsc0.0";
				}
			}
			if (getSettings()
					.isWriteDomainArchitectureSimilarityScoresToOutput()) {
				csvRow += buildDomainWeightColumn(bestBr);
			}
		}
		return csvRow;
	}
}