package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastSearchContentAdapter;

public class BlastSearchContentAdapterTest {

	private BlastSearchContentAdapter blastSearchContentAdapter;

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
		this.blastSearchContentAdapter = new BlastSearchContentAdapter(
				TestUtils.mockProteinDb(), "swissprot");
	}

	@Test
	public void testPassesBlacklist() {
		List<String> blacklist = getSettings().getBlastResultsBlackList(
				"swissprot");
		blacklist.add("^Probable");
		blacklist.add("^Putative");
		blacklist.add("similar");
		String blastResultDescriptionLine = "Probable sheep hormone";
		assertTrue(!this.blastSearchContentAdapter
				.passesBlacklist(blastResultDescriptionLine));
		blastResultDescriptionLine = "Putative goat cheese";
		assertTrue(!this.blastSearchContentAdapter
				.passesBlacklist(blastResultDescriptionLine));
		blastResultDescriptionLine = "Protein similar to goats and sheep.";
		assertTrue(!this.blastSearchContentAdapter
				.passesBlacklist(blastResultDescriptionLine));
		blastResultDescriptionLine = "This is a truly important Description-Line a bio-informatician wants to pass the blacklist.";
		assertTrue(this.blastSearchContentAdapter
				.passesBlacklist(blastResultDescriptionLine));
	}

	@Test
	public void testFilter() {
		List<String> filter = getSettings().getBlastResultsFilter("swissprot");
		filter.add("\\sOS=\\S+");
		filter.add("\\shappy\\snew\\syear");
		filter.add("bat\\d{3}sheep\\s");
		String blastResultDescriptionLine = "Probable sheep hormone OS=Baaaa";
		assertEquals("Probable sheep hormone", this.blastSearchContentAdapter
				.filter(blastResultDescriptionLine));

		blastResultDescriptionLine = "Probable sheep hormone happy new year";
		assertEquals("Probable sheep hormone", this.blastSearchContentAdapter
				.filter(blastResultDescriptionLine));

		blastResultDescriptionLine = "Probable bat666sheep hormone";
		assertEquals("Probable hormone", this.blastSearchContentAdapter
				.filter(blastResultDescriptionLine));
	}
}
