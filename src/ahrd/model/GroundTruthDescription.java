package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import java.util.HashSet;
import java.util.Set;

public class GroundTruthDescription {

	private Set<String> tokens = new HashSet<String>();
	private String accession;
	private String description;

	public GroundTruthDescription() {
		super();
	}

	public static GroundTruthDescription constructFromFastaEntry(String fastaEntry) {
		GroundTruthDescription rd = new GroundTruthDescription();

		// First line is a combination of Accession and Description
		String[] fastaData = fastaEntry.split("\n");
		// First token before whitespace-char is the Accession
		rd.setAccession(fastaData[0].split(" ")[0].trim());
		// Everything after the Accession is considered the description-line:
		rd.setDescription(fastaData[0].replace(rd.getAccession(), "").trim());
		// Process the ground truth's human readable description as requested by
		// the user (Settings) -
		// NOTE, if the HRD passes the Blacklist and no filtering is
		// requested the HRD does not have to be processed any further.
		if (getSettings().getGroundTruthDescriptionBlacklist() != null
				&& !getSettings().getGroundTruthDescriptionBlacklist().isEmpty()) {
			if (!DescriptionScoreCalculator.passesBlacklist(rd.getDescription(),
					getSettings().getGroundTruthDescriptionBlacklist())) {
				// Does NOT pass blacklist
				rd.setDescription("");
			} else if (getSettings().getGroundTruthDescriptionFilter() != null
					&& !getSettings().getGroundTruthDescriptionFilter().isEmpty()) {
				// Passes Blacklist AND is requested to be filtered:
				rd.setDescription(DescriptionScoreCalculator.filter(rd.getDescription(),
						getSettings().getGroundTruthDescriptionFilter()));
			}
		}
		// Tokenize, and if requested in Settings retain only those tokens that
		// pass the Blacklist:
		rd.setTokens(TokenScoreCalculator.tokenize(rd.getDescription(), getSettings().getGroundTruthTokenBlacklist()));
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
