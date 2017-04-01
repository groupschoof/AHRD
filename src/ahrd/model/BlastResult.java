package ahrd.model;

import static ahrd.controller.Settings.getSettings;

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

import ahrd.controller.Settings;
import ahrd.exception.MissingProteinException;

/**
 * This class merges attributes from BlastHits: accession, description and
 * Blast's High Scoring Pair: start, end, bitScore, eValue
 * 
 * @author klee, hallab
 */
public class BlastResult implements Comparable<BlastResult> {

	public static final String TOKEN_SPLITTER_REGEX = "-|/|;|\\\\|,|:|\"|'|\\.|\\s+|\\||\\(|\\)";
	public static final String FASTA_PROTEIN_HEADER_ACCESSION_GROUP_NAME = "accession";
	public static final String FASTA_PROTEIN_HEADER_DESCRIPTION_GROUP_NAME = "description";
	public static final String SHORT_ACCESSION_GROUP_NAME = "shortAccession";
	public static final String GO_TERM_GROUP_NAME = "goTerm";

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
	private Fscore evaluationScore;
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
	 * Gene Ontology term annotations. Populated in the best blast results of
	 * proteins if the output of the best blast results and the evaluation of GO
	 * terms is requested.
	 * Or populated if the output of the highest possible GO annotation score is requested. 
	 */
	private Set<GOterm> goAnnotations = new HashSet<GOterm>();
	/**
	 * GO annotation F scores
	 */
	private Fscore simpleGoAnnotationScore; // Simple cardinality based
	private Fscore ancestryGoAnnotationScore; // Cardinality of the term ancestries 
	private Fscore semSimGoAnnotationScore; // Based on semantic similarity derived from term information content 

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
	 * @param uniqueShortAccessions
	 *            - Used only if AHRD is requested to generate Gene Ontology
	 *            term annotations
	 * @throws IOException
	 * @throws MissingProteinException
	 */
	public static void readBlastResults(Map<String, Protein> proteinDb, String blastDbName,
			Set<String> uniqueAccessions) throws MissingProteinException, IOException {
		Map<String, List<BlastResult>> brs = parseBlastResults(proteinDb, blastDbName, uniqueAccessions);
		parseBlastDatabase(proteinDb, blastDbName, brs);
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
	 * @param uniqueShortAccessions
	 *            - Used only if AHRD is requested to generate Gene Ontology
	 *            term annotations
	 * @return Map<String,List<BlastResult>> Set of Hit-Accessions (Key) to the
	 *         full BlastResult(s) (Value)
	 * @throws MissingProteinException
	 * @throws IOException
	 */
	public static Map<String, List<BlastResult>> parseBlastResults(Map<String, Protein> proteinDb, String blastDbName,
			Set<String> uniqueShortAccessions) throws MissingProteinException, IOException {
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
					addBlastResult(brs, br, uniqueShortAccessions);
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
	 * @param uniqueShortAccessions
	 */
	public static void addBlastResult(Map<String, List<BlastResult>> brs, BlastResult br,
			Set<String> uniqueShortAccessions) {
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
		// Finally, if AHRD is requested to annotate Gene Ontology Terms, we
		// need to extract all unique short reference protein (BlastResult)
		// accessions:
		if (getSettings().hasGeneOntologyAnnotations()) {
			uniqueShortAccessions.add(br.getShortAccession());
		}
	}

	/**
	 * Adds the sequence length and Human Readable Description (HRD) to all
	 * matching BlastHits found in the argument Map 'blastResults'. Afterwards
	 * the respective BlastResult instances are added to their respective
	 * Proteins. See function <code>generateHRDCandidateForProtein</code> for
	 * more details.
	 * 
	 * @param blastResults
	 * @param fastaAccession
	 * @param hitAALength
	 * @param hrd
	 */
	public static void fastaEntryValuesForBlastHit(Map<String, List<BlastResult>> blastResults, String fastaAccession,
			Integer hitAALength, String hrd) {
		for (BlastResult br : blastResults.get(fastaAccession)) {
			br.setSubjectLength(hitAALength);
			br.setDescription(hrd);
			br.generateHRDCandidateForProtein();
		}
	}

	/**
	 * Extracts from the provided protein database in FASTA format the length
	 * and human readable descriptions of those Proteins that are Hits
	 * (Subjects) in the argument blastResults. Each time such a Hit is found
	 * the mentioned measurements are set in the respective instance of
	 * BlastResult and subsequently the method 'Protein.addBlastResult' is
	 * invoked.
	 * 
	 * @param proteinDb
	 * @param blastDbName
	 * @param blastResults
	 * @throws IOException
	 */
	public static void parseBlastDatabase(Map<String, Protein> proteinDb, String blastDbName,
			Map<String, List<BlastResult>> blastResults) throws IOException {
		// Parse line by line FASTA Blast search DB. Extract Subject Lengths and
		// Subject HRDs.
		BufferedReader fastaIn = null;
		try {
			fastaIn = new BufferedReader(new FileReader(getSettings().getPathToBlastDatabase(blastDbName)));
			String str, hrd = new String();
			String acc = "";
			Integer hitAALength = new Integer(0);
			boolean hit = false;
			while ((str = fastaIn.readLine()) != null) {
				if (str.startsWith(">")) {
					// Finished reading in the original Fasta-Entry of a
					// Blast-Hit? If so, process it:
					if (hit) {
						fastaEntryValuesForBlastHit(blastResults, acc, hitAALength, hrd);
						// Clean up to enable processing the next Hit
						hitAALength = new Integer(0);
						// Note, that the boolean 'hit' will be set in the
						// following If-Else-Block.

					}

					// Process the current Fasta-Header-Line:
					Matcher m = getSettings().getFastaHeaderRegex(blastDbName).matcher(str);
					if (!m.matches()) {
						// Provided REGEX to parse FASTA header does not work in
						// this case:
						System.err.println("WARNING: FASTA header line\n" + str.trim()
								+ "\ndoes not match provided regular expression\n"
								+ getSettings().getFastaHeaderRegex(blastDbName).toString()
								+ "\n. The header and the following entry, including possibly respective matching BLAST Hits, are ignored and discarded.\n"
								+ "To fix this, please use - Blast database specific - parameter "
								+ Settings.FASTA_HEADER_REGEX_KEY
								+ " to provide a regular expression that matches ALL FASTA headers in Blast database '"
								+ blastDbName + "'.");
					} else if (blastResults.containsKey(m.group(FASTA_PROTEIN_HEADER_ACCESSION_GROUP_NAME).trim())) {
						// Found the next Blast HIT:
						acc = m.group(FASTA_PROTEIN_HEADER_ACCESSION_GROUP_NAME).trim();
						hrd = m.group(FASTA_PROTEIN_HEADER_DESCRIPTION_GROUP_NAME).trim();
						// Following lines, until the next header, contain
						// information to be collected:
						hit = true;
					} else {
						// Found a Protein in the FASTA database, that is of no
						// relevance within this context:
						hit = false;
					}
				} else if (hit) {
					// Process non header-line, if and only if, we are reading
					// the sequence of a Blast-Hit:
					hitAALength += str.trim().length();
				}
			}
			// Was the last read FASTA entry a Blast-Hit? If so, it needs
			// processing:
			if (hit)
				fastaEntryValuesForBlastHit(blastResults, acc, hitAALength, hrd);
		} finally {
			fastaIn.close();
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
		Set<String> tknBlackList = getSettings().getTokenBlacklist(getBlastDatabaseName());
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
		Set<String> blacklist = getSettings().getBlastResultsBlacklist(getBlastDatabaseName());
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
	 * is obtained from. You can evaluate AHRD based <i>only</i> on tokens that
	 * passed the Blacklist and Filtering with the correct input parameter. See
	 * <code>Settings.evaluateValidTokens</code> for details.
	 * 
	 * @Note: This method uses the static respective static method from
	 *        Model-Class GroundTruthDescription.
	 */
	public void tokenizeForEvaluation() {
		if (getSettings().getEvaluateOnlyValidTokens())
			setEvaluationTokens(getTokens());
		else
			setEvaluationTokens(TokenScoreCalculator.tokenize(getDescription(), new HashSet<String>()));
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
			theClone.setTokens(TokenScoreCalculator.tokenize(theClone.getDescription(), new HashSet<String>()));
			getProtein().getEvaluationScoreCalculator().addBestUnchangedBlastResult(getBlastDatabaseName(), theClone);
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

	public Fscore getEvaluationScore() {
		return evaluationScore;
	}

	public void setEvaluationScore(Fscore evaluationScore) {
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

	public Set<GOterm> getGoAnnotations() {
		return goAnnotations;
	}

	public void setGoAnnotations(Set<GOterm> goAnnotation) {
		this.goAnnotations = goAnnotation;
	}

	public Fscore getSimpleGoAnnotationScore() {
		return simpleGoAnnotationScore;
	}

	public void setSimpleGoAnnotationScore(Fscore simpleGoAnnotationScore) {
		this.simpleGoAnnotationScore = simpleGoAnnotationScore;
	}

	public Fscore getAncestryGoAnnotationScore() {
		return ancestryGoAnnotationScore;
	}

	public void setAncestryGoAnnotationScore(Fscore ancestryGoAnnotationScore) {
		this.ancestryGoAnnotationScore = ancestryGoAnnotationScore;
	}

	public Fscore getSemSimGoAnnotationScore() {
		return semSimGoAnnotationScore;
	}

	public void setSemSimGoAnnotationScore(Fscore semSimGoAnnotationScore) {
		this.semSimGoAnnotationScore = semSimGoAnnotationScore;
	}

}
