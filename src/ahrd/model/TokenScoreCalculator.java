package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.roundToNDecimalPlaces;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Functions calculating Token-Scores.
 */
public class TokenScoreCalculator {

	private Map<String, Double> cumulativeTokenBitScores = new HashMap<String, Double>();
	private Map<String, Double> cumulativeTokenBlastDatabaseScores = new HashMap<String, Double>();
	private Map<String, Double> cumulativeTokenOverlapScores = new HashMap<String, Double>();
	private double totalTokenBitScore = 0;
	private double totalTokenBlastDatabaseScore = 0;
	private double totalTokenOverlapScore = 0;
	private Map<String, Double> tokenScores = new HashMap<String, Double>();
	private Protein protein;
	// Please enter your initials ___
	private double tokenHighScore = 0.0;

	public static boolean tokenPassesBlacklist(String token, List<String> blacklist) {
		// No Token passes being NULL or empty String
		boolean passed = (token != null && !token.equals(""));
		Iterator<String> i = blacklist.iterator();
		while (i.hasNext() && passed) {
			Pattern p = Pattern.compile(i.next());
			Matcher m = p.matcher(token);
			// A Match against a blacklisted RegExp lets the token fail:
			passed = !m.find();
		}
		return passed;
	}

	/**
	 * Splits the String description into tokens using the regular expressions
	 * found in <code>BlastResult.TOKEN _SPLITTER_REGEX</code> and then passes
	 * them through Blacklisting retaining only tokens that pass.
	 * 
	 * @param description
	 * @param tokenBlacklist
	 * @return An instance of Set holding all tokens extracted from description
	 */
	public static Set<String> tokenize(String description, List<String> tokenBlacklist) {
		Set<String> tkns = new HashSet<String>();
		for (String tokenCandidate : new HashSet<String>(
				Arrays.asList(description.split(BlastResult.TOKEN_SPLITTER_REGEX)))) {
			tokenCandidate = tokenCandidate.toLowerCase();
			if (tokenCandidate != null && !tokenCandidate.equals("")
					&& tokenPassesBlacklist(tokenCandidate, tokenBlacklist))
				tkns.add(tokenCandidate);
		}
		return tkns;
	}

	/**
	 * Computes the quality of the local sequence alignment between the Query
	 * and a Subject (see BLAST for details). This quality is assessed as the
	 * proportion of the sum of aligned residues divided by the cumulative
	 * sequence length.
	 * 
	 * @param queryStart
	 * @param queryEnd
	 * @param queryLength
	 * @param subjectStart
	 * @param subjectEnd
	 * @param subjectLength
	 * @return double
	 */
	public static double overlapScore(double queryStart, double queryEnd, double queryLength, double subjectStart,
			double subjectEnd, double subjectLength) {
		return ((queryEnd - queryStart + 1.0) + (subjectEnd - subjectStart + 1.0)) / (queryLength + subjectLength);
	}

	public TokenScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	/**
	 * Returns the sum of BlastResult's Tokens' Scores.
	 */
	public double descriptionLineSummedTokenScore(BlastResult br) {
		double sum = 0.0;
		for (String token : br.getTokens()) {
			sum += getTokenScores().get(token);
		}
		return sum;
	}

	/**
	 * Assigns each Token in each BlastResult it's TokenScore and stores it in
	 * the Map 'tokenScores' key is Token and value is it's TokenScore.
	 */
	public void assignTokenScores() {
		// iterate through blast databases having BlastResults
		for (String iterBlastDb : getProtein().getBlastResults().keySet()) {
			// iterate through blast results in coming from the different blast
			// databases
			for (BlastResult iterResult : getProtein().getBlastResults().get(iterBlastDb)) {
				// iterate through tokens in different blast result desc-lines
				for (String token : iterResult.getTokens()) {
					if (!(getTokenScores().containsKey(token))) {
						double tokenscore = tokenScore(token, iterBlastDb);
						getTokenScores().put(token, new Double(tokenscore));
						// remember highest token score
						if (tokenscore > getTokenHighScore()) {
							setTokenHighScore(tokenscore);
						}
					}
				}
			}
		}
	}

	/**
	 * Iterates over all Tokens in the Map 'tokenScores' and re-assigns them new
	 * TokenScores. Each non-informative Token is assigned the new TokenScore :=
	 * (old TokenScore) - (tokenHighScore / 2).
	 */
	public void filterTokenScores() {
		for (String token : getTokenScores().keySet()) {
			if (!isInformativeToken(token)) {
				getTokenScores().put(token, new Double(getTokenScores().get(token) - getTokenHighScore() / 2));
			}
		}
	}

	/**
	 * Informative tokens have a token-score greater than half of the
	 * tokenHighScore.
	 * 
	 * @Note: This method cannot be invoked before having the tokenScores-Set
	 *        initialized!
	 */
	public boolean isInformativeToken(String token) {
		return getTokenScores().get(token) > getTokenHighScore() / 2;
	}

	/**
	 * Once per BlastResult's unique token the following <em>cumulative</em>
	 * scores are measured: 1. BitScore 2. DatabaseScore 3. OverlapScore
	 * 
	 * @param BlastResult
	 *            br
	 */
	public void measureCumulativeScores(BlastResult br) {
		for (String token : br.getTokens()) {
			Double overlapScore = TokenScoreCalculator.overlapScore(br.getQueryStart(), br.getQueryEnd(),
					getProtein().getSequenceLength(), br.getSubjectStart(), br.getSubjectEnd(), br.getSubjectLength());
			addCumulativeTokenBitScore(token, br.getBitScore());
			addCumulativeTokenBlastDatabaseScore(token, br.getBlastDatabaseName());
			addCumulativeTokenOverlapScore(token, overlapScore);
		}
	}

	/**
	 * Once per BlastResult the following <em>total</em> scores are measured: 1.
	 * BitScore 2. DatabaseScore 3. OverlapScore
	 * 
	 * @param BlastResult
	 *            br
	 */
	public void measureTotalScores(BlastResult br) {
		Double overlapScore = TokenScoreCalculator.overlapScore(br.getQueryStart(), br.getQueryEnd(),
				getProtein().getSequenceLength(), br.getSubjectStart(), br.getSubjectEnd(), br.getSubjectLength());
		setTotalTokenBlastDatabaseScore(
				getTotalTokenBlastDatabaseScore() + getSettings().getBlastDbWeight(br.getBlastDatabaseName()));
		setTotalTokenOverlapScore(getTotalTokenOverlapScore() + overlapScore);
		setTotalTokenBitScore(getTotalTokenBitScore() + br.getBitScore());
	}

	/**
	 * @param token
	 * @return token-score
	 */
	public double tokenScore(String token, String blastDatabaseName) {
		// Validate:
		Double bitScoreWeight = getSettings().getTokenScoreBitScoreWeight();
		Double databaseScoreWeight = getSettings().getTokenScoreDatabaseScoreWeight();
		Double overlapScoreWeight = getSettings().getTokenScoreOverlapScoreWeight();
		double validateSumToOne = roundToNDecimalPlaces(bitScoreWeight + databaseScoreWeight + overlapScoreWeight, 9);
		// Tolerate rounding error <= 10^-3
		if (!(validateSumToOne >= 0.999 && validateSumToOne <= 1.001))
			throw new IllegalArgumentException(
					"The three weights 'bitScoreWeight', 'databaseScoreWeight', and 'overlapScoreWeight' should sum up to 1, but actually sum up to: "
							+ (bitScoreWeight + databaseScoreWeight + overlapScoreWeight));
		// Calculate Token-Score:
		return (bitScoreWeight * getCumulativeTokenBitScores().get(token) / getTotalTokenBitScore()
				+ databaseScoreWeight * getCumulativeTokenBlastDatabaseScores().get(token)
						/ getTotalTokenBlastDatabaseScore()
				+ overlapScoreWeight * getCumulativeTokenOverlapScores().get(token) / getTotalTokenOverlapScore());
	}

	public void addCumulativeTokenBitScore(String token, double bitScore) {
		if (!getCumulativeTokenBitScores().containsKey(token))
			getCumulativeTokenBitScores().put(token, new Double(bitScore));
		else
			getCumulativeTokenBitScores().put(token, new Double(bitScore + getCumulativeTokenBitScores().get(token)));
	}

	public void addCumulativeTokenOverlapScore(String token, double overlapScore) {
		if (!getCumulativeTokenOverlapScores().containsKey(token))
			getCumulativeTokenOverlapScores().put(token, new Double(overlapScore));
		else
			getCumulativeTokenOverlapScores().put(token,
					new Double(overlapScore + getCumulativeTokenOverlapScores().get(token)));
	}

	public void addCumulativeTokenBlastDatabaseScore(String token, String blastDatabaseName) {
		Integer blastDatabaseWeight = getSettings().getBlastDbWeight(blastDatabaseName);
		if (!getCumulativeTokenBlastDatabaseScores().containsKey(token))
			getCumulativeTokenBlastDatabaseScores().put(token, new Double(blastDatabaseWeight));
		else
			getCumulativeTokenBlastDatabaseScores().put(token,
					new Double(blastDatabaseWeight + getCumulativeTokenBlastDatabaseScores().get(token)));
	}

	public double sumOfAllTokenScores(BlastResult blastResult) {
		double sum = 0.0;
		for (String token : blastResult.getTokens()) {
			sum += getTokenScores().get(token);
		}
		return sum;
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

	public Map<String, Double> getCumulativeTokenBitScores() {
		return cumulativeTokenBitScores;
	}

	public void setCumulativeTokenBitScores(Map<String, Double> cumulativeTokenBitScores) {
		this.cumulativeTokenBitScores = cumulativeTokenBitScores;
	}

	public Map<String, Double> getCumulativeTokenBlastDatabaseScores() {
		return cumulativeTokenBlastDatabaseScores;
	}

	public void setCumulativeTokenBlastDatabaseScores(Map<String, Double> cumulativeTokenBlastDatabaseScores) {
		this.cumulativeTokenBlastDatabaseScores = cumulativeTokenBlastDatabaseScores;
	}

	public Map<String, Double> getCumulativeTokenOverlapScores() {
		return cumulativeTokenOverlapScores;
	}

	public void setCumulativeTokenOverlapScores(Map<String, Double> cumulativeTokenOverlapScores) {
		this.cumulativeTokenOverlapScores = cumulativeTokenOverlapScores;
	}

	public double getTotalTokenBitScore() {
		return totalTokenBitScore;
	}

	public void setTotalTokenBitScore(double totalTokenBitScore) {
		this.totalTokenBitScore = totalTokenBitScore;
	}

	public double getTotalTokenBlastDatabaseScore() {
		return totalTokenBlastDatabaseScore;
	}

	public void setTotalTokenBlastDatabaseScore(double totalTokenBlastDatabaseScore) {
		this.totalTokenBlastDatabaseScore = totalTokenBlastDatabaseScore;
	}

	public double getTotalTokenOverlapScore() {
		return totalTokenOverlapScore;
	}

	public void setTotalTokenOverlapScore(double totalTokenOverlapScore) {
		this.totalTokenOverlapScore = totalTokenOverlapScore;
	}

	/**
	 * Get tokenHighScore.
	 * 
	 * @return tokenHighScore as double.
	 */
	public double getTokenHighScore() {
		return tokenHighScore;
	}

	/**
	 * Set tokenHighScore.
	 * 
	 * @param tokenHighScore
	 *            the value to set.
	 */
	public void setTokenHighScore(double tokenHighScore) {
		this.tokenHighScore = tokenHighScore;
	}

	/**
	 * Get tokenScores.
	 * 
	 * @return tokenScores as Map<String, Double>.
	 */
	public Map<String, Double> getTokenScores() {
		return tokenScores;
	}

	/**
	 * Set tokenScores.
	 * 
	 * @param tokenScores
	 *            the value to set.
	 */
	public void setTokenScores(Map<String, Double> tokenScores) {
		this.tokenScores = tokenScores;
	}
}
