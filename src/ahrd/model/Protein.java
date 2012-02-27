package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ahrd.exception.MissingAccessionException;

public class Protein {
	
	private String accession;
	private Integer sequenceLength;
	private Map<String, List<BlastResult>> blastResults;
	private Set<InterproResult> interproResults = new HashSet<InterproResult>();
	private Set<GeneOntologyResult> goResults = new HashSet<GeneOntologyResult>();
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
		// set to train AHRD's parameters, a EvaluationScoreCalculator is needed:
		if (getSettings().getWriteBestBlastHitsToOutput()
				|| getSettings().isInTrainingMode())
			setEvaluationScoreCalculator(new EvaluationScoreCalculator(this));
	}

	public static Protein constructFromFastaEntry(String fastaEntry)
			throws MissingAccessionException {
		String[] fasta_data = fastaEntry.split("\n");
		String accession = fasta_data[0].split(" ")[0];
		if (accession == null || accession.equals("")) {
			throw new MissingAccessionException(
					"Missing protein-accession in:\n" + fastaEntry);
		}
		String[] sequence_parts = new String[fasta_data.length - 1];
		System.arraycopy(fasta_data, 1, sequence_parts, 0,
				fasta_data.length - 1);
		String sequence = "";
		for (String sequence_part : sequence_parts) {
			sequence += sequence_part.trim();
		}
		return new Protein(accession, sequence.length());
	}

	/**
	 * Construct Memory-Database of Proteins!
	 * 
	 * @param fastaFileContent
	 * @return
	 */
	public static Map<String, Protein> initializeProteins(
			String fastaFileContent) throws MissingAccessionException {
		Map<String, Protein> proteins = new HashMap<String, Protein>();
		String[] fastaEntries = fastaFileContent.trim().split(">");
		for (String fastaEntry : fastaEntries) {
			if (fastaEntry != null && !fastaEntry.trim().equals("")
					&& !fastaEntry.equals("\n")) {
				Protein prot = constructFromFastaEntry(fastaEntry);
				proteins.put(prot.accession, prot);
			}
		}
		return proteins;
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
		// Measure the frequency of the BlastResult-Description-Line:
		getDescriptionScoreCalculator().measureDescriptionLineFrequency(
				br.patternize());
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

	public Set<GeneOntologyResult> getGoResults() {
		return goResults;
	}

	public void setGoResults(Set<GeneOntologyResult> goResults) {
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
}
