package ahrd.model;

import static ahrd.controller.Settings.ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.DESCRIPTION_GROUP_NAME;
import static ahrd.controller.Settings.DEFAULT_LINE_SEP;
import static ahrd.controller.Settings.getSettings;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;

import ahrd.exception.MissingAccessionException;

public class GroundTruthDescription {

	private Set<String> tokens = new HashSet<String>();
	private String accession;
	private String description;

	public GroundTruthDescription() {
		super();
	}

	public static GroundTruthDescription constructFromFastaEntry(String fastaEntry) throws MissingAccessionException {
		GroundTruthDescription rd = new GroundTruthDescription();
		
		String[] fasta_data = fastaEntry.split(DEFAULT_LINE_SEP);
		Matcher m = getSettings().getGroundTruthFastaRegex().matcher(fasta_data[0]);
		if (m.find()) {
			rd.setAccession(m.group(ACCESSION_GROUP_NAME));
			if (rd.getAccession() == null || rd.getAccession().equals("")) {
				throw new MissingAccessionException("Missing protein-accession in:\n" + fastaEntry);
			} else {
				rd.setDescription(m.group(DESCRIPTION_GROUP_NAME));
			}
		} else {
			throw new MissingAccessionException("Missing protein-accession in:\n" + fastaEntry);
		}
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
