package ahrd.model;

import static ahrd.controller.Settings.DEFAULT_LINE_SEP;
import static ahrd.controller.Settings.ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.getSettings;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import ahrd.exception.MissingAccessionException;

public class Protein {

	private String accession;
	private String sequence;
	private Integer sequenceLength;
	private Map<String, List<BlastResult>> blastResults;
	private Set<InterproResult> interproResults = new HashSet<InterproResult>();
	private Set<String> goResults = new HashSet<String>();
	private Set<GOterm> goResultsTerms = new HashSet<GOterm>();
	private Set<GOterm> goSlimTerms = new HashSet<GOterm>();
	private Map<GOterm, Double> goResultsTermsConfidence = new HashMap<GOterm, Double>();
	private Map<String, Double> goCentricTermConfidences = new HashMap<String, Double>();
	private TokenScoreCalculator tokenScoreCalculator;
	private LexicalScoreCalculator lexicalScoreCalculator;
	private DescriptionScoreCalculator descriptionScoreCalculator;
	private EvaluationScoreCalculator evaluationScoreCalculator;

	public Protein(String accession, Integer sequenceLength) {
		super();
		setAccession(accession);
		setSequenceLength(sequenceLength);
		setBlastResults(new HashMap<String, List<BlastResult>>());
		setTokenScoreCalculator(new TokenScoreCalculator(this));
		setLexicalScoreCalculator(new LexicalScoreCalculator(this));
		setDescriptionScoreCalculator(new DescriptionScoreCalculator(this));
		// If current Settings require evaluation of AHRD's performance or are
		// set to train AHRD's parameters, a EvaluationScoreCalculator is
		// needed:
		if (getSettings().getWriteBestBlastHitsToOutput()
				|| getSettings().isInTrainingMode())
			setEvaluationScoreCalculator(new EvaluationScoreCalculator(this));
	}

	public Protein(String accession, String aaSequence) {
		super();
		setAccession(accession);
		setSequence(aaSequence);
		setSequenceLength(aaSequence.length());
		setBlastResults(new HashMap<String, List<BlastResult>>());
		setTokenScoreCalculator(new TokenScoreCalculator(this));
		setLexicalScoreCalculator(new LexicalScoreCalculator(this));
		setDescriptionScoreCalculator(new DescriptionScoreCalculator(this));
		// If current Settings require evaluation of AHRD's performance or are
		// set to train AHRD's parameters, a EvaluationScoreCalculator is
		// needed:
		if (getSettings().getWriteBestBlastHitsToOutput()
				|| getSettings().isInTrainingMode())
			setEvaluationScoreCalculator(new EvaluationScoreCalculator(this));
	}

	public static List<String> splitFasta(String fastaStr) {
		List<String> fastaEntries = new ArrayList<String>(Arrays.asList(fastaStr.split("(^|\n|\r)(?=>)")));
		fastaEntries.removeAll(Arrays.asList("", null));
		return fastaEntries;
	}

	public static Protein constructFromFastaEntry(String fastaEntry) throws MissingAccessionException {
		String[] fasta_data = fastaEntry.split(DEFAULT_LINE_SEP);
		Matcher m = getSettings().getProteinsFastaRegex().matcher(fasta_data[0]);
		String acc = "";
		if (m.find()) {
			acc = m.group(ACCESSION_GROUP_NAME);
			if (acc == null || acc.equals("")) {
				throw new MissingAccessionException("Missing protein-accession in:\n" + fastaEntry);
			}
		} else {
			throw new MissingAccessionException("Missing protein-accession in:\n" + fastaEntry);
		}
		String[] sequence_parts = new String[fasta_data.length - 1];
		System.arraycopy(fasta_data, 1, sequence_parts, 0,
				fasta_data.length - 1);
		String sequence = "";
		for (String sequence_part : sequence_parts) {
			sequence += sequence_part.trim();
		}
		// Construct the new Protein, either storing its AA-sequence or just the
		// sequence's length:
		Protein p = null;
		if (getSettings().doOutputFasta())
			p = new Protein(acc, sequence);
		else
			p = new Protein(acc, sequence.length());
		return p;
	}

	/**
	 * Construct Memory-Database of Proteins!
	 * 
	 * @param fastaFileContent
	 * @return
	 */
	public static Map<String, Protein> initializeProteins(String fastaFileContent) throws MissingAccessionException {
		Map<String, Protein> proteins = new HashMap<String, Protein>();
		List<String> fastaEntries = splitFasta(fastaFileContent);
		for (String fastaEntry : fastaEntries) {
			if (fastaEntry != null && !fastaEntry.trim().equals("")	&& !fastaEntry.equals("\n")) {
				Protein prot = constructFromFastaEntry(fastaEntry);
				proteins.put(prot.accession, prot);
			}
		}
		return proteins;
	}

	/**
	 * Extracts all unique Gene Ontology (GO) terms annotated to the Proteins in
	 * argument prots.
	 * 
	 * @param prots
	 * @return Set<String>
	 */
	public static Set<String> uniqueGOaccessions(Collection<Protein> prots) {
		Set<String> ugt = new HashSet<String>();
		for (Protein prot : prots) {
			ugt.addAll(prot.getGoResults());
		}
		return ugt;
	}

	/**
	 * Adds the BlastResult to the Protein's set and measures the cumulative and
	 * total scores later needed to calculate the Token-Scores. Also finds the
	 * highest BitScore and Description-Line-Frequency. The argument BlastResult
	 * is expected to have passed Blacklist and Filter and is expected to have
	 * been token- and patternized.
	 * 
	 * @param BlastResult
	 */
	public void addBlastResult(BlastResult br) {
		String blastDb = br.getBlastDatabaseName();
		if (!getBlastResults().containsKey(blastDb)) {
			getBlastResults().put(blastDb, new ArrayList<BlastResult>());
		}
		getBlastResults().get(blastDb).add(br);
		// Measure TokenScore related cumulative Scores:
		getTokenScoreCalculator().measureCumulativeScores(br);
		// Measure TokenScore related total Scores:
		getTokenScoreCalculator().measureTotalScores(br);
		// Measure highest BitScore:
		getDescriptionScoreCalculator().measureMaxBitScore(br.getBitScore());
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public Map<String, List<BlastResult>> getBlastResults() {
		return blastResults;
	}

	public void setBlastResults(Map<String, List<BlastResult>> blastResults) {
		this.blastResults = blastResults;
	}

	public Set<InterproResult> getInterproResults() {
		return interproResults;
	}

	public void setInterproResults(Set<InterproResult> interproResults) {
		this.interproResults = interproResults;
	}

	public Set<String> getGoResults() {
		return goResults;
	}

	public void setGoResults(Set<String> goResults) {
		this.goResults = goResults;
	}

	public TokenScoreCalculator getTokenScoreCalculator() {
		return tokenScoreCalculator;
	}

	public void setTokenScoreCalculator(
			TokenScoreCalculator tokenScoreCalculator) {
		this.tokenScoreCalculator = tokenScoreCalculator;
	}

	/**
	 * Get lexicalScoreCalculator.
	 * 
	 * @return lexicalScoreCalculator as LexicalScoreCalculator.
	 */
	public LexicalScoreCalculator getLexicalScoreCalculator() {
		return lexicalScoreCalculator;
	}

	/**
	 * Set lexicalScoreCalculator.
	 * 
	 * @param lexicalScoreCalculator
	 *            the value to set.
	 */
	public void setLexicalScoreCalculator(
			LexicalScoreCalculator lexicalScoreCalculator) {
		this.lexicalScoreCalculator = lexicalScoreCalculator;
	}

	/**
	 * Get descriptionScoreCalculator.
	 * 
	 * @return descriptionScoreCalculator as DescriptionScoreCalculator.
	 */
	public DescriptionScoreCalculator getDescriptionScoreCalculator() {
		return descriptionScoreCalculator;
	}

	/**
	 * Set descriptionScoreCalculator.
	 * 
	 * @param descriptionScoreCalculator
	 *            the value to set.
	 */
	public void setDescriptionScoreCalculator(
			DescriptionScoreCalculator descriptionScoreCalculator) {
		this.descriptionScoreCalculator = descriptionScoreCalculator;
	}

	/**
	 * Get sequenceLength.
	 * 
	 * @return sequenceLength as Integer.
	 */
	public Integer getSequenceLength() {
		return sequenceLength;
	}

	/**
	 * Set sequenceLength.
	 * 
	 * @param sequenceLength
	 *            the value to set.
	 */
	public void setSequenceLength(Integer sequenceLength) {
		this.sequenceLength = sequenceLength;
	}

	public EvaluationScoreCalculator getEvaluationScoreCalculator() {
		return evaluationScoreCalculator;
	}

	public void setEvaluationScoreCalculator(
			EvaluationScoreCalculator evaluationScoreCalculator) {
		this.evaluationScoreCalculator = evaluationScoreCalculator;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public Set<GOterm> getGoResultsTerms() {
		return goResultsTerms;
	}

	public void setGoResultsTerms(Set<GOterm> goResultsTerms) {
		this.goResultsTerms = goResultsTerms;
	}

	public Set<GOterm> getGoSlimTerms() {
		return goSlimTerms;
	}

	public void setGoSlimTerms(Set<GOterm> goSlimTerms) {
		this.goSlimTerms = goSlimTerms;
	}

	public Map<GOterm, Double> getGoResultsTermsConfidence() {
		return goResultsTermsConfidence;
	}

	public void setGoResultsTermsConfidence(Map<GOterm, Double> goResultsTermsConfidence) {
		this.goResultsTermsConfidence = goResultsTermsConfidence;
	}

	public Map<String, Double> getGoCentricTermConfidences() {
		return goCentricTermConfidences;
	}

	public void setGoCentricTermConfidences(Map<String, Double> goCentricTermConfidences) {
		this.goCentricTermConfidences = goCentricTermConfidences;
	}

}
