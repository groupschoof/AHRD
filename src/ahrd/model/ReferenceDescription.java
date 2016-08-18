package ahrd.model;

import static ahrd.controller.Settings.getSettings;
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
		// Process the reference's human readable description as requested by
		// the user (Settings) -
		// NOTE, if the HRD passes the Blacklist and no filtering is
		// requested the HRD does not have to be processed any further.
		if (getSettings().getReferencesDescriptionBlacklist() != null
				&& !getSettings().getReferencesDescriptionBlacklist().isEmpty()) {
			if (!DescriptionScoreCalculator.passesBlacklist(rd.getDescription(),
					getSettings().getReferencesDescriptionBlacklist())) {
				// Does NOT pass blacklist
				rd.setDescription("");
			} else if (getSettings().getReferencesDescriptionFilter() != null
					&& !getSettings().getReferencesDescriptionFilter().isEmpty()) {
				// Passes Blacklist AND is requested to be filtered:
				rd.setDescription(DescriptionScoreCalculator.filter(rd.getDescription(),
						getSettings().getReferencesDescriptionFilter()));
			}
		}
		// Tokenize, and if requested in Settings retain only those tokens that
		// pass the Blacklist:
		rd.setTokens(TokenScoreCalculator.tokenize(rd.getDescription(), getSettings().getReferencesTokenBlacklist()));
		return rd;
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
