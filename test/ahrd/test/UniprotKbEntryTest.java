package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.HashSet;

import nu.xom.ParsingException;
import nu.xom.ValidityException;

import org.junit.Test;

import ahrd.model.UniprotKBEntry;

public class UniprotKbEntryTest {

	private static final String ACCESSION = "Q0KFR8";

	@Test
	public void testUrl() throws UnsupportedEncodingException {
		assertEquals(
				"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/Q0KFR8/xml",
				UniprotKBEntry.url(ACCESSION));
		assertEquals(
				"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/sp%7CP34529%7CDCR1_CAEEL/xml",
				UniprotKBEntry.url("sp|P34529|DCR1_CAEEL"));
	}

	@Test
	public void testFromUrl() throws IOException, ValidityException,
			ParsingException {
		UniprotKBEntry u = UniprotKBEntry
				.fromUrl(UniprotKBEntry.url(ACCESSION));
		assertNotNull("Expected an instance of UniprotKBEntry but got NULL.", u);
		assertEquals(ACCESSION, u.getAccession());
		assertEquals(
				new HashSet<String>(Arrays.asList(new String[] { "IPR020591",
						"IPR013159", "IPR024633", "IPR010921", "IPR003593",
						"IPR018312", "IPR001957", "IPR013317" })),
				u.getIprAnnotations());
		assertEquals(
				new HashSet<String>(Arrays.asList(new String[] { "PF00308",
						"PF08299", "PF11638" })), u.getPfamAnnotations());
	}
}
