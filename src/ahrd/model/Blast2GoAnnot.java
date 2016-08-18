package ahrd.model;

import static ahrd.model.TokenScoreCalculator.tokenize;

import java.util.ArrayList;
import java.util.Set;

public class Blast2GoAnnot implements Comparable<Blast2GoAnnot> {

	private String accession;
	private String description;
	private Set<String> evaluationTokens;
	private Double evaluationScore = 0.0;

	public static Blast2GoAnnot fromBlast2GoEntry(String resultLine) {
		Blast2GoAnnot res = null;
		String[] vals = resultLine.split("\t");
		String accession = vals[0].trim();
		// GO-Term-Accession is in position 2, which is ignored here.
		String description = vals[2].trim();
		if (accession != null && description != null && !accession.equals("") && !description.equals(""))
			res = new Blast2GoAnnot(accession, description);
		return res;
	}

	public Blast2GoAnnot(String accession, String description) {
		super();
		setAccession(accession);
		setDescription(description);
		setEvaluationTokens(tokenize(getDescription(), new ArrayList<String>()));
	}

	/**
	 * Compares to argument Blast2GoAnnot by comparing the appropriate
	 * Evaluation-Scores.
	 */
	public int compareTo(Blast2GoAnnot b2ga) {
		return this.getEvaluationScore().compareTo(b2ga.getEvaluationScore());
	}

	/**
	 * Blast2GoAnnots shall equal whenever their descriptions is the same.
	 * 
	 * @Note: Remember that this is only safe to use, when in the context of a
	 *        single and unique <em>accession</em>.
	 */
	@Override
	public int hashCode() {
		return getDescription().hashCode();
	}

	/**
	 * Blast2GoAnnots shall equal whenever their descriptions is the same.
	 * 
	 * @Note: Remember that this is only safe to use, when in the context of a
	 *        single and unique <em>accession</em>.
	 */
	@Override
	public boolean equals(Object o) {
		// check for self-comparison
		if (this == o)
			return true;
		// object of different class is not equal to this
		if (!(o instanceof Blast2GoAnnot))
			return false;
		// Do the respective description equal one another?
		return ((Blast2GoAnnot) o).getDescription().equals(getDescription());
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public Set<String> getEvaluationTokens() {
		return evaluationTokens;
	}

	public void setEvaluationTokens(Set<String> evaluationTokens) {
		this.evaluationTokens = evaluationTokens;
	}

	public Double getEvaluationScore() {
		return evaluationScore;
	}

	public void setEvaluationScore(Double evaluationScore) {
		this.evaluationScore = evaluationScore;
	}
}
