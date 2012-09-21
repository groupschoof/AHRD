package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.model.ReferenceDescription.tokenizeDescription;
import static ahrd.model.TokenScoreCalculator.passesBlacklist;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.program.sax.BlastLikeSAXParser;
import org.biojava.bio.program.ssbind.SeqSimilarityAdapter;
import org.biojava.bio.search.SearchContentHandler;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

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
	private List<Double> domainWeights;
	/**
	 * Based on the algorithm proposed in
	 * "Protein comparison at the domain architecture level" by Lee and Lee
	 * (http://www.biomedcentral.com/1471-2105/10/S15/S5) we compute the
	 * BlastResult's domain similarity score as the similarity of its domain
	 * architecture to the query protein's.
	 * 
	 * @Note: The score ranges from 0.0 to 1.0 and is the higher the more
	 *        similar the two domain architectures are.
	 */
	private Double domainSimilarityScore;

	public BlastResult(String blastDatabaseName) {
		super();
		setBlastDatabaseName(blastDatabaseName);
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
	 * Reads in BlastResults and assigns them to the Proteins in argument
	 * proteinDb.
	 * 
	 * Example code taken from biojava.org.
	 */
	public static void parseBlastResults(Map<String, Protein> proteinDb,
			String blastDbName) throws MissingProteinException, SAXException,
			IOException {
		// get the Blast input as a Stream
		InputStream is = new FileInputStream(getSettings()
				.getPathToBlastResults(blastDbName));
		// make a BlastLikeSAXParser
		BlastLikeSAXParser parser = new BlastLikeSAXParser();
		// try to parse, even if the blast version is not recognized.
		parser.setModeLazy();
		// make the SAX event adapter that will pass events to a Handler.
		SeqSimilarityAdapter adapter = new SeqSimilarityAdapter();
		// set the parsers SAX event adapter
		parser.setContentHandler(adapter);
		// register builder with custom adapter
		SearchContentHandler scHandler = new BlastSearchContentAdapter(
				proteinDb, blastDbName);
		adapter.setSearchContentHandler(scHandler);
		// parse the file, after this the result List will be populated with
		// SeqSimilaritySearchResults
		parser.parse(new InputSource(is));
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
		for (String tokenCandidate : new HashSet<String>(Arrays
				.asList(getDescription().split(TOKEN_SPLITTER_REGEX)))) {
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

	public List<Double> getDomainWeights() {
		return domainWeights;
	}

	public void setDomainWeights(List<Double> domainWeights) {
		this.domainWeights = domainWeights;
	}

	public Double getDomainSimilarityScore() {
		return domainSimilarityScore;
	}

	public void setDomainSimilarityScore(Double domainSimilarityScore) {
		this.domainSimilarityScore = domainSimilarityScore;
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
}
