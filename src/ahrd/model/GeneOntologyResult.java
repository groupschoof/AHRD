package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import java.util.Map;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;

public class GeneOntologyResult implements Comparable<GeneOntologyResult> {

	private String acc;
	private String name;
	private double probability;

	public GeneOntologyResult(String acc, String name, double probability) {
		super();
		setAcc(acc);
		setName(name);
		setProbability(probability);
	}

	public static void parseGeneOntologyResult(Map<String, Protein> proteinDb)
			throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(
				getSettings().getPathToGeneOntologyResults())));
		String iterLine = null;
		while ((iterLine = br.readLine()) != null) {
			String[] parts = iterLine.split("\t");
			Protein prot = proteinDb.get(parts[0].trim());
			if (prot != null) {
				double goProbability = Double.parseDouble((parts[1]));
				String goAcc = parts[2].trim();
				if (goAcc != null && !goAcc.equals("")) {
					GeneOntologyResult gr = new GeneOntologyResult(goAcc,
							parts[3].trim(), goProbability);
					prot.getGoResults().add(gr);
				}
			}
		}
	}

	public int compareTo(GeneOntologyResult goToCompare) {
		return this.getAcc().compareTo(goToCompare.getAcc());
	}

	/**
	 * Get acc.
	 * 
	 * @return acc as String.
	 */
	public String getAcc() {
		return acc;
	}

	/**
	 * Set acc.
	 * 
	 * @param acc
	 *            the value to set.
	 */
	public void setAcc(String acc) {
		this.acc = acc;
	}

	/**
	 * Get name.
	 * 
	 * @return name as String.
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set name.
	 * 
	 * @param name
	 *            the value to set.
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Get probability.
	 * 
	 * @return probability as double.
	 */
	public double getProbability() {
		return probability;
	}

	/**
	 * Set probability.
	 * 
	 * @param probability
	 *            the value to set.
	 */
	public void setProbability(double probability) {
		this.probability = probability;
	}
}
