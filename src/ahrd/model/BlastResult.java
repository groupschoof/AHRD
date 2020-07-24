package ahrd.model;

import static ahrd.controller.Settings.SHORT_ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.getSettings;
import static ahrd.model.AhrdDb.getReferenceProteinDAO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;

/**
 * This class merges attributes from BlastHits: accession, description and
 * Blast's High Scoring Pair: start, end, bitScore, eValue
 * 
 * @author klee, hallab
 */
public class BlastResult implements Comparable<BlastResult> {

	private String accession;
	private String shortAccession;
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
	 * Makes a double string representation parseable by Double.parseDouble
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
	public BlastResult(String accession, double eValue, int queryStart, int queryEnd, int subjectStart, int subjectEnd,
			double bitScore, String blastDatabaseName, Protein protein) {
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

	public BlastResult(String accession, double eValue, String description, int queryStart, int queryEnd,
			int subjectStart, int subjectEnd, int subjectLength, double bitScore, String blastDatabaseName) {
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
	 * Wraps up two steps:
	 * <ul>
	 * <li>Parse tabular Sequence Similarity Search Results</li>
	 * <li>Extract Human Readable Descriptions (HRDs) and Sequence Lengths from
	 * Protein Database in FASTA format</li>
	 * </ul>
	 * 
	 * @param proteinDb
	 * @param blastDbName
	 * @throws IOException
	 * @throws MissingProteinException
	 * @throws MissingAccessionException
	 */
	public static void readBlastResults(Map<String, Protein> proteinDb, String blastDbName)
			throws MissingProteinException, IOException, MissingAccessionException {
		Map<String, List<BlastResult>> brs = parseBlastResults(proteinDb, blastDbName);
		parseBlastDatabase(brs);
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
	public static Map<String, List<BlastResult>> parseBlastResults(Map<String, Protein> proteinDb, String blastDbName)
			throws MissingProteinException, IOException {
		Map<String, List<BlastResult>> brs = new HashMap<String, List<BlastResult>>();
		BufferedReader fastaIn = null;
		try {
			fastaIn = new BufferedReader(new FileReader(getSettings().getPathToBlastResults(blastDbName)));
			String str;
			while ((str = fastaIn.readLine()) != null) {
				// Only evaluate current line, either if there is no
				// comment-line-regex given, or if it is given AND it does not
				// match:
				if (getSettings().getSeqSimSearchTableCommentLineRegex() == null
						|| !getSettings().getSeqSimSearchTableCommentLineRegex().matcher(str).matches()) {
					String[] brFields = str.split(getSettings().getSeqSimSearchTableSep());
					if (!proteinDb.containsKey(brFields[getSettings().getSeqSimSearchTableQueryCol()])) {
						throw new MissingProteinException("Could not find Protein for Accession '"
								+ brFields[getSettings().getSeqSimSearchTableQueryCol()] + "' in Protein Database.");
					} // ELSE
					BlastResult br = new BlastResult(brFields[getSettings().getSeqSimSearchTableSubjectCol()],
							Double.parseDouble(validateDouble(brFields[getSettings().getSeqSimSearchTableEValueCol()])),
							Integer.parseInt(brFields[getSettings().getSeqSimSearchTableQueryStartCol()]),
							Integer.parseInt(brFields[getSettings().getSeqSimSearchTableQueryEndCol()]),
							Integer.parseInt(brFields[getSettings().getSeqSimSearchTableSubjectStartCol()]),
							Integer.parseInt(brFields[getSettings().getSeqSimSearchTableSubjectEndCol()]),
							Double.parseDouble(brFields[getSettings().getSeqSimSearchTableBitScoreCol()]), blastDbName,
							proteinDb.get(brFields[getSettings().getSeqSimSearchTableQueryCol()]));
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
	public static void addBlastResult(Map<String, List<BlastResult>> brs, BlastResult br) {
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

	/**
	 * Extracts from AHRD's persistent database the length and human readable
	 * descriptions of those Proteins that are Hits (Subjects) in the argument
	 * blastResults. Each time such a Hit is found the mentioned measurements
	 * are set in the respective instance of BlastResult and subsequently the
	 * method 'Protein.addBlastResult' is invoked.
	 * 
	 * @param blastResults
	 * @throws MissingAccessionException
	 */
	public static void parseBlastDatabase(Map<String, List<BlastResult>> blastResults)
			throws MissingAccessionException {
		ReferenceProtein rp;
		for (String accession : blastResults.keySet()) {
			rp = getReferenceProteinDAO().byAccession.get(accession);
			if (rp == null) {
				throw new MissingAccessionException("Found Blast-Hits to reference protein '" + accession
						+ "' but could not find the matching Reference-Protein in AHRD's persistent Database."
						+ " Please update it accordingly.");
			}
			setReferenceProteinValuesInBlastHits(rp, blastResults.get(accession));
		}
	}

	/**
	 * Populates all BlastResults (Hits) referring to the argument reference
	 * protein <code>rp</code> with information extracted from the original
	 * sequence database in Fasta-Format, i.e. the description and the sequence
	 * (Hit) total length. After doing so prepares the Hit to be processed by
	 * AHRD, i.e. passes it through the Blacklist, Filtering, and Tokenizing
	 * (see <code>generateHRDCandidateForProtein</code> for more details).
	 * 
	 * @param rp
	 * @param hits
	 */
	public static void setReferenceProteinValuesInBlastHits(ReferenceProtein rp, List<BlastResult> hits) {
		for (BlastResult blastResult : hits) {
			blastResult.setDescription(rp.getHrd());
			blastResult.setSubjectLength(rp.getSequenceLength());
			blastResult.generateHRDCandidateForProtein();
		}
	}

	public static List<BlastResult> filterBestScoringBlastResults(List<BlastResult> blastResults, int howMany) {
		if (blastResults.size() > howMany) {
			List<BlastResult> sortedBlastResults = new ArrayList<BlastResult>(blastResults);
			Collections.sort(sortedBlastResults);
			blastResults = sortedBlastResults.subList(0, howMany);
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

	/**
	 * Splits the Human Readable Description into single tokens and stores them
	 * in this' tokens field.
	 */
	public void tokenize() {
		List<String> tknBlackList = getSettings().getTokenBlackList(getBlastDatabaseName());
		this.setTokens(TokenScoreCalculator.tokenize(this.getDescription(), tknBlackList));
	}

	/**
	 * Checks if this' human readable description passes the blacklist
	 * associated with this' sequence database name.
	 * 
	 * @param blastResultDescriptionLine
	 * @return boolean TRUE if and only if the human readable description passes
	 *         the respective blacklist. FALSE otherwise.
	 */
	public boolean passesBlacklist(String blastResultDescriptionLine) {
		List<String> blacklist = getSettings().getBlastResultsBlackList(getBlastDatabaseName());
		return DescriptionScoreCalculator.passesBlacklist(blastResultDescriptionLine, blacklist);
	}

	/**
	 * Filters this' human readable description using the global filter
	 * implemented in <code>DescriptionScoreCalculator.filter(...)</code>
	 * 
	 * @param blastResultDescriptionLine
	 * @return String the modified version of this' human readable description,
	 *         in which all matches to the respective filters are deleted.
	 */
	public String filter(String blastResultDescriptionLine) {
		List<String> filter = getSettings().getBlastResultsFilter(getBlastDatabaseName());
		return DescriptionScoreCalculator.filter(blastResultDescriptionLine, filter);
	}

	/**
	 * Evaluation Tokens are <i>not</i> filtered with the TOKEN-BLACKLIST, as we
	 * want to evaluate <i>all</i> tokens, that are printed out, too. This set
	 * of evaluation-tokens is set if and only if, AHRD is run in Evaluator-Mode
	 * and this BlastResult is the best scoring of the Blast-Search-Result, it
	 * is obtained from. You can evaluate AHRD based <i>only<\i> on tokens that
	 * passed the Blacklist and Filtering with the correct input parameter. See
	 * <code>Settings.evaluateValidTokens</code> for details.
	 * 
	 * @Note: This method uses the static respective static method from
	 *        Model-Class ReferenceDescription.
	 */
	public void tokenizeForEvaluation() {
		if (getSettings().getEvaluateValidTokens())
			setEvaluationTokens(getTokens());
		else
			setEvaluationTokens(TokenScoreCalculator.tokenize(getDescription(), new ArrayList<String>()));
	}

	public boolean isValid() {
		return (getAccession() != null && (!getAccession().equals("")) && getBitScore() != null
				&& getDescription() != null && (!getDescription().equals("")) && getQueryEnd() != null
				&& getQueryStart() != null && (getQueryStart() < getQueryEnd()) && getSubjectEnd() != null
				&& getSubjectStart() != null && (getSubjectEnd() > getSubjectStart()) && getSubjectLength() != null
				&& getEValue() != null && getTokens() != null && getTokens().size() > 0
				&& getBlastDatabaseName() != null
				&& getSettings().getBlastDatabases().contains(getBlastDatabaseName()));
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
		return new BlastResult(new String(this.getAccession()), new Double(eValue), new String(description),
				new Integer(queryStart), new Integer(queryEnd), new Integer(subjectStart), new Integer(subjectEnd),
				new Integer(subjectLength), new Double(bitScore), new String(blastDatabaseName));
	}

	/**
	 * Investigates this instance's properties, especially the Description. If
	 * the instance is valid and its description passes the Blacklist, it will
	 * be added as a candidate HRD to the respective query Protein's
	 * BlastResults.
	 */
	public void generateHRDCandidateForProtein() {
		// For Training-Purposes:
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			// Of course we do have to treat this best-blast-hit
			// differently than the further to process one below, so
			// clone:
			BlastResult theClone = clone();
			// Pass best Blast-Hit's Description through filter:
			theClone.setDescription(filter(theClone.getDescription()));
			// Tokenize without filtering tokens through the Blacklist:
			theClone.setTokens(TokenScoreCalculator.tokenize(theClone.getDescription(), new ArrayList<String>()));
			getProtein().getEvaluationScoreCalculator().addUnchangedBlastResult(getBlastDatabaseName(), theClone);
		}
		if (passesBlacklist(getDescription())) {
			// Pass bestScoringHSP through filter:
			setDescription(filter(getDescription()));
			// Tokenize the filtered Description-Line:
			tokenize();
			// Pass bestScoringHSP through Blacklist and add it, if it
			// is still valid:
			if (isValid()) {
				// Generate the pattern of the Description's unique
				// tokens:
				patternize();
				// Adds the BlastResult to the getProtein()'s set and
				// measures the cumulative and total scores later needed
				// to calculate the Token-Scores:
				getProtein().addBlastResult(this);
			}
		}
	}

	/**
	 * Extracts from the possibly longer Accession the shorter one which is used
	 * in the reference Gene Ontology Annotation (GOA) file. This is required
	 * due to the fact that UniprotKB uses short accessions in the GOA files but
	 * provides long accessions in the Blast Databases. A very unfortunate lack
	 * of standardization, indeed!
	 * 
	 * @return String
	 */
	public String getShortAccession() {
		if (shortAccession == null) {
			Pattern p = getSettings().getShortAccessionRegex(getBlastDatabaseName());
			Matcher m = p.matcher(getAccession());
			setShortAccession(getAccession());
			if (!m.find()) {
				System.err.println("WARNING: Regular Expression '" + p.toString()
						+ "' does NOT match - using pattern.find(...) - Blast Hit Accession '" + getAccession()
						+ "' - continuing with the original accession. This might lead to unrecognized reference GO annotations!");
			} else {
				shortAccession = m.group(SHORT_ACCESSION_GROUP_NAME);
			}
		}
		return (shortAccession);
	}

	public void setShortAccession(String shortAccession) {
		this.shortAccession = shortAccession;
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
			throw new IllegalArgumentException("Blast-Database-Name must not be NULL.");
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
