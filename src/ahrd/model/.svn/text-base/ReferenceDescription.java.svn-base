package ahrd.model;

import java.util.HashSet;
import java.util.Set;

public class ReferenceDescription {

	private Set<String> tokens = new HashSet<String>();
	private String accession;
	private String description;

	public ReferenceDescription() {
		super();
	}

	public static ReferenceDescription constructFromFastaEntry(String fastaEntry) {
		ReferenceDescription rd = new ReferenceDescription();

		// First line is a combination of Accession and Description
		String[] fastaData = fastaEntry.split("\n");
		// First token before whitespace-char is the Accession
		rd.setAccession(fastaData[0].split(" ")[0].trim());
		// Everything after the Accession is considered the description-line:
		rd.setDescription(fastaData[0].replace(rd.getAccession(), "").trim());
		// Tokenize description and filter out tokens matching any regex in the
		// blacklist:
		rd.setTokens(tokenizeDescription(rd.getDescription()));

		return rd;
	}

	/**
	 * Tokenizes a String using Blastresult.TOKEN_SPLITTER_REGEX and returns all
	 * resulting unique tokens.
	 * 
	 * @param description
	 * @param blacklist
	 * @return Set<String>
	 */
	public static Set<String> tokenizeDescription(String description) {
		Set<String> tkns = new HashSet<String>();
		for (String tkn : description.split(BlastResult.TOKEN_SPLITTER_REGEX)) {
			String tokenCandidate = tkn.trim().toLowerCase();
			if (tokenCandidate != null && !tokenCandidate.equals(""))
				tkns.add(tokenCandidate);
		}
		return tkns;
	}

	public Set<String> getTokens() {
		return tokens;
	}

	public void setTokens(Set<String> tokens) {
		this.tokens = tokens;
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

}
