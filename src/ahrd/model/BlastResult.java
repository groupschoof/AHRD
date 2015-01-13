package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.model.ReferenceDescription.tokenizeDescription;
import static ahrd.model.TokenScoreCalculator.passesBlacklist;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ahrd.exception.MissingProteinException;

/**
 * This class merges attributes from BlastHits: accession, description and
 * Blast's High Scoring Pair: start, end, bitScore, eValue
 * 
 * @author klee, hallab
 */
public class BlastResult implements Comparable<BlastResult> {

	public static final String TOKEN_SPLITTER_REGEX = "-|/|;|\\\\|,|:|\"|'|\\.|\\s+|\\||\\(|\\)";

	private String accession;
	private Double eValue;
	private String description;
	/**
	 * Query's start position in local alignment
	 */
	private Integer queryStart;
	/**
	 * Query's stop position in local alignment
	 */
	private Integer queryEnd;
	/**
	 * Subject's start position in local alignment
	 */
	private Integer subjectStart;
	/**
	 * Subject's stop position in local alignment
	 */
	private Integer subjectEnd;
	/**
	 * Length of subject's amino acid sequence.
	 */
	private Integer subjectLength;
	private Double bitScore;
	private String blastDatabaseName;
	/**
	 * The descriptionScore is calculated by AHRD.
	 */
	private Double descriptionScore;
	private Set<String> tokens = new HashSet<String>();
	/**
	 * The evaluationScore is calculated while training or evaluating AHRD's
	 * performance in comparison with the "Best Blast Hit"-Method:
	 */
	private Double evaluationScore;
	/**
	 * Evaluation Tokens are <i>not</i> filtered with the TOKEN-BLACKLIST, as we
	 * want to evaluate <i>all</i> tokens, that are printed out, too. This set
	 * of evaluation-tokens is set if and only if, AHRD is run in Evaluator-Mode
	 * and this BlastResult is the best scoring of the Blast-Search-Result, it
	 * is obtained from.
	 */
	private Set<String> evaluationTokens;
	/**
	 * The query accession is only stored during the parsing of tabular sequence
	 * similarity search results. It should only be used in that context.
	 */
	private Protein protein;

	/**
	 * Makes a double string representation parseable by Double.parseDouble(…).
	 * 
	 * @param input
	 * @return String
	 */
	public static String validateDouble(String input) {
		if (input.startsWith("e") || input.startsWith("E"))
			input = "1" + input;
		return input;
	}

	public BlastResult(String blastDatabaseName) {
		super();
		setBlastDatabaseName(blastDatabaseName);
	}

	/**
	 * Constructor missing only Subject-Length and the Human Readable
	 * Description, both of which will be obtained from the original Database in
	 * FASTA format.
	 * 
	 * @param accession
	 * @param eValue
	 * @param queryStart
	 * @param queryEnd
	 * @param subjectStart
	 * @param subjectEnd
	 * @param bitScore
	 * @param blastDatabaseName
	 * @param protein
	 * 
	 * @return BlastResult
	 */
	public BlastResult(String accession, double eValue, int queryStart,
			int queryEnd, int subjectStart, int subjectEnd, double bitScore,
			String blastDatabaseName, Protein protein) {
		super();
		setAccession(accession);
		setEValue(eValue);
		setQueryStart(queryStart);
		setQueryEnd(queryEnd);
		setSubjectStart(subjectStart);
		setSubjectEnd(subjectEnd);
		setBitScore(bitScore);
		setBlastDatabaseName(blastDatabaseName);
		setProtein(protein);
	}

	public BlastResult(String accession, double eValue, String description,
			int queryStart, int queryEnd, int subjectStart, int subjectEnd,
			int subjectLength, double bitScore, String blastDatabaseName) {
		super();
		setAccession(accession);
		setEValue(eValue);
		setDescription(description);
		setQueryStart(queryStart);
		setQueryEnd(queryEnd);
		setSubjectStart(subjectStart);
		setSubjectEnd(subjectEnd);
		setSubjectLength(subjectLength);
		setBitScore(bitScore);
		setBlastDatabaseName(blastDatabaseName);
	}

	public BlastResult(String blastDbName, String accession, String description) {
		super();
		setBlastDatabaseName(blastDbName);
		setAccession(accession);
		setDescription(description);
	}

	/**
	 * Reads in BlastResults, or any other results from sequence similarity
	 * searches, and assigns them to the Proteins in argument proteinDb. The
	 * result file is expected to be in tabular format, that is each line is
	 * supposed to contain a single High Scoring Pair (HSP). Preferred format is
	 * 'Blast8' (-m 8).
	 * 
	 * @param proteinDb
	 * @param blastDbName
	 * @return Map<String,List<BlastResult>> Set of Hit-Accessions (Key) to the
	 *         full BlastResult(s) (Value)
	 * @throws MissingProteinException
	 * @throws IOException
	 */
	public static Map<String, List<BlastResult>> parseBlastResults(
			Map<String, Protein> proteinDb, String blastDbName)
			throws MissingProteinException, IOException {
		Map<String, List<BlastResult>> brs = new HashMap<String, List<BlastResult>>();
		BufferedReader fastaIn = null;
		try {
			fastaIn = new BufferedReader(new FileReader(getSettings()
					.getPathToBlastResults(blastDbName)));
			String str;
			while ((str = fastaIn.readLine()) != null) {
				// Only evaluate current line, either if there is no
				// comment-line-regex given, or if it is given AND it does not
				// match:
				if (getSettings().getSeqSimSearchTableCommentLineRegex() == null
						|| !getSettings()
								.getSeqSimSearchTableCommentLineRegex()
								.matcher(str).matches()) {
					String[] brFields = str.split(getSettings()
							.getSeqSimSearchTableSep());
					if (!proteinDb.containsKey(brFields[getSettings()
							.getSeqSimSearchTableQueryCol()])) {
						throw new MissingProteinException(
								"Could not find Protein for Accession '"
										+ brFields[getSettings()
												.getSeqSimSearchTableQueryCol()]
										+ "' in Protein Database.");
					}// ELSE
					BlastResult br = new BlastResult(
							brFields[getSettings()
									.getSeqSimSearchTableSubjectCol()],
							Double.parseDouble(validateDouble(brFields[getSettings()
									.getSeqSimSearchTableEValueCol()])),
							Integer.parseInt(brFields[getSettings()
									.getSeqSimSearchTableQueryStartCol()]),
							Integer.parseInt(brFields[getSettings()
									.getSeqSimSearchTableQueryEndCol()]),
							Integer.parseInt(brFields[getSettings()
									.getSeqSimSearchTableSubjectStartCol()]),
							Integer.parseInt(brFields[getSettings()
									.getSeqSimSearchTableSubjectEndCol()]),
							Double.parseDouble(brFields[getSettings()
									.getSeqSimSearchTableBitScoreCol()]),
							blastDbName, proteinDb.get(brFields[getSettings()
									.getSeqSimSearchTableQueryCol()]));
					addBlastResult(brs, br);
				}
			}
		} finally {
			fastaIn.close();
		}
		return brs;
	}

	/**
	 * The argument BlastResult is added to the argument Map of BlastResults. If
	 * a BlastResult of same accession and for the same query protein is already
	 * present, and the argument BlastResult has a better Bit-Score, it replaces
	 * the one of worse Bit-Score. If this is a so far not seen BlastResult,
	 * regarding Hit-Accession and Query-Accession, it will simply be added.
	 * 
	 * @param brs
	 * @param br
	 */
	public static void addBlastResult(Map<String, List<BlastResult>> brs,
			BlastResult br) {
		if (brs.containsKey(br.getAccession())) {
			boolean isMultipleHsp = false;
			List<BlastResult> sameHitBrs = brs.get(br.getAccession());
			for (BlastResult iterBr : sameHitBrs) {
				// If and only if there is a BlastResult of same Hit and Query
				// that also has a better Bit-Score, replace the one of lower
				// Bit-Score with the higher one:
				if (iterBr.getProtein().equals(br.getProtein())) {
					isMultipleHsp = true;
					if (iterBr.getBitScore() < br.getBitScore()) {
						sameHitBrs.remove(iterBr);
						sameHitBrs.add(br);
					}
				}
			}
			// If this a Hit for another Protein, add it:
			if (!isMultipleHsp) {
				sameHitBrs.add(br);
			}

		} else {
			// Add a new List<BlastResult> containing the argument BlastResult
			// br to the argument Map brs:
			List<BlastResult> sameHitBrs = new ArrayList<BlastResult>();
			sameHitBrs.add(br);
			brs.put(br.getAccession(), sameHitBrs);
		}
	}

	public static void parseBlastDatabase(Map<String, Protein> proteinDb,
			String blastDbName, Map<String, List<BlastResult>> blastResults)
			throws IOException {
		// Parse line by line FASTA Blast search DB. Extract Subject Lengths and
		// Subject HRDs.
		BufferedReader fastaIn = null;
		try {
			fastaIn = new BufferedReader(new FileReader(getSettings()
					.getPathToBlastDatabase(blastDbName)));
			String str;
			while ((str = fastaIn.readLine()) != null) {
				if (str.startsWith(">")) {

				}
			}
		} finally {
			fastaIn.close();
		}
	}

	public static List<BlastResult> filterBestScoringBlastResults(
			List<BlastResult> blastResults, int howMany) {
		if (blastResults.size() > howMany) {
			List<BlastResult> sortedBlastResults = new ArrayList<BlastResult>(
					blastResults);
			Collections.sort(sortedBlastResults);
			blastResults = sortedBlastResults.subList(sortedBlastResults.size()
					- howMany - 1, sortedBlastResults.size() - 1);
		}
		return blastResults;
	}

	/**
	 * Sorts unique tokens and returns them concatenated. @NOTE: As this is
	 * expectedly fast, we do not need to store the generated pattern in the
	 * BlastResult. We can just generate it again, when needed.
	 * 
	 * @return String pattern
	 */
	public String patternize() {
		String pattern = "";
		if (getTokens().size() > 0) {
			List<String> sortedTkns = new ArrayList<String>(getTokens());
			Collections.sort(sortedTkns);
			for (String tkn : sortedTkns) {
				pattern += tkn;
			}
		}
		return pattern;
	}

	public void tokenize() {
		List<String> tknBlackList = getSettings().getTokenBlackList(
				getBlastDatabaseName());
		for (String tokenCandidate : new HashSet<String>(
				Arrays.asList(getDescription().split(TOKEN_SPLITTER_REGEX)))) {
			tokenCandidate = tokenCandidate.toLowerCase();
			if (passesBlacklist(tokenCandidate, tknBlackList))
				getTokens().add(tokenCandidate);
		}
	}

	/**
	 * Evaluation Tokens are <i>not</i> filtered with the TOKEN-BLACKLIST, as ee
	 * want to evaluate <i>all</i> tokens, that are printed out, too. This set
	 * of evaluation-tokens is set if and only if, AHRD is run in Evaluator-Mode
	 * and this BlastResult is the best scoring of the Blast-Search-Result, it
	 * is obtained from.
	 * 
	 * @Note: This method uses the static respective static method from
	 *        Model-Class ReferenceDescription.
	 */
	public void tokenizeForEvaluation() {
		setEvaluationTokens(tokenizeDescription(getDescription()));
	}

	public boolean isValid() {
		return (getAccession() != null && (!getAccession().equals(""))
				&& getBitScore() != null && getDescription() != null
				&& (!getDescription().equals("")) && getQueryEnd() != null
				&& getQueryStart() != null && (getQueryStart() < getQueryEnd())
				&& getSubjectEnd() != null && getSubjectStart() != null
				&& (getSubjectEnd() > getSubjectStart())
				&& getSubjectLength() != null && getEValue() != null
				&& getTokens() != null && getTokens().size() > 0 && getSettings()
				.getBlastDatabases().contains(getBlastDatabaseName()));
	}

	/**
	 * Compares to argument BlastResult by comparing the appropriate E-Values.
	 */
	public int compareTo(BlastResult compareBlastResult) {
		return this.getEValue().compareTo(compareBlastResult.getEValue());
	}

	/**
	 * Dolly was a sweet little sheep!
	 * 
	 * @return A clone of this Instance, cloning all attributes except
	 *         descriptionScore, tokens and evaluationScore.
	 */
	public BlastResult clone() {
		return new BlastResult(new String(this.getAccession()), new Double(
				eValue), new String(description), new Integer(queryStart),
				new Integer(queryEnd), new Integer(subjectStart), new Integer(
						subjectEnd), new Integer(subjectLength), new Double(
						bitScore), new String(blastDatabaseName));
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public Double getEValue() {
		return eValue;
	}

	public void setEValue(Double value) {
		eValue = value;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public Integer getQueryStart() {
		return queryStart;
	}

	public void setQueryStart(Integer start) {
		this.queryStart = start;
	}

	public Integer getQueryEnd() {
		return queryEnd;
	}

	public void setQueryEnd(Integer end) {
		this.queryEnd = end;
	}

	public Double getBitScore() {
		return bitScore;
	}

	public void setBitScore(Double bitScore) {
		this.bitScore = bitScore;
	}

	public Set<String> getTokens() {
		return tokens;
	}

	public void setTokens(Set<String> tokens) {
		this.tokens = tokens;
	}

	public Double getDescriptionScore() {
		return descriptionScore;
	}

	public void setDescriptionScore(Double descriptionScore) {
		this.descriptionScore = descriptionScore;
	}

	public String getBlastDatabaseName() {
		return blastDatabaseName;
	}

	public void setBlastDatabaseName(String blastDatabaseName) {
		if (blastDatabaseName == null)
			throw new IllegalArgumentException(
					"Blast-Database-Name must not be NULL.");
		this.blastDatabaseName = blastDatabaseName;
	}

	public Double getEvaluationScore() {
		return evaluationScore;
	}

	public void setEvaluationScore(Double evaluationScore) {
		this.evaluationScore = evaluationScore;
	}

	public Set<String> getEvaluationTokens() {
		return evaluationTokens;
	}

	public void setEvaluationTokens(Set<String> evaluationTokens) {
		this.evaluationTokens = evaluationTokens;
	}

	public Integer getSubjectStart() {
		return subjectStart;
	}

	public void setSubjectStart(Integer subjectStart) {
		this.subjectStart = subjectStart;
	}

	public Integer getSubjectEnd() {
		return subjectEnd;
	}

	public void setSubjectEnd(Integer subjectEnd) {
		this.subjectEnd = subjectEnd;
	}

	public Integer getSubjectLength() {
		return subjectLength;
	}

	public void setSubjectLength(Integer subjectLength) {
		this.subjectLength = subjectLength;
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

}
