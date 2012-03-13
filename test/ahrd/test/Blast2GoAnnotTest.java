package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import ahrd.model.Blast2GoAnnot;

public class Blast2GoAnnotTest {

	@Test
	public void testFromBlast2GoEntry() {
		String b2gResEntry = "accession\tGO:123456\tSheep horn growthase";
		Blast2GoAnnot b2ga = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		// test:
		assertEquals("accession", b2ga.getAccession());
		assertEquals("Sheep horn growthase", b2ga.getDescription());
		assertTrue(
				"Expected Token in Blast2GoAnnots Token-Set, but did not find it.",
				b2ga.getEvaluationTokens().contains("sheep"));
		assertTrue(
				"Expected Token in Blast2GoAnnots Token-Set, but did not find it.",
				b2ga.getEvaluationTokens().contains("horn"));
		assertTrue(
				"Expected Token in Blast2GoAnnots Token-Set, but did not find it.",
				b2ga.getEvaluationTokens().contains("growthase"));
	}

	@Test
	public void testEquals() {
		String b2gResEntry = "accession\tGO:123456\tSheep horn growthase";
		String b2gResEntryTwo = "accession\tGO:654321\tGoat wool grothase";
		Blast2GoAnnot b2gaRef = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		Blast2GoAnnot b2gaClone = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		Blast2GoAnnot b2gaTwo = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntryTwo);
		// test
		assertTrue("Expecting Blast2GoAnnot to equal itself.", b2gaRef
				.equals(b2gaRef));
		assertTrue(
				"Expecting Blast2GoAnnot to equal one with the same description.",
				b2gaRef.equals(b2gaClone));
		assertTrue(
				"Expecting Blast2GoAnnot NOT to equal Object of different class.",
				!b2gaRef.equals("object of different class"));
		assertTrue("Expecting Blast2GoAnnot NOT to equal null.", !b2gaRef
				.equals(null));
		assertTrue(
				"Expecting Blast2GoAnnot NOT to equal Blast2GoAnnot of different description.",
				!b2gaRef.equals(b2gaTwo));
	}

	@Test
	public void testCompare() {
		String b2gResEntry = "accession\tGO:123456\tSheep horn growthase";
		String b2gResEntryTwo = "accession\tGO:654321\tGoat wool grothase";
		Blast2GoAnnot b2gaRef = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		b2gaRef.setEvaluationScore(0.1);
		Blast2GoAnnot b2gaClone = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		b2gaClone.setEvaluationScore(0.1);
		Blast2GoAnnot b2gaTwo = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntryTwo);
		b2gaTwo.setEvaluationScore(0.9);
		// test
		assertEquals(-1, b2gaRef.compareTo(b2gaTwo));
		assertEquals(1, b2gaTwo.compareTo(b2gaRef));
		assertEquals(0, b2gaRef.compareTo(b2gaClone));
		assertEquals(0, b2gaRef.compareTo(b2gaRef));
	}

	@Test
	public void testHashCode() {
		String b2gResEntry = "accession\tGO:123456\tSheep horn growthase";
		String b2gResEntryTwo = "accession\tGO:654321\tGoat wool grothase";
		Blast2GoAnnot b2gaRef = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		Blast2GoAnnot b2gaClone = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		Blast2GoAnnot b2gaTwo = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntryTwo);
		// test
		assertTrue(
				"Expecting hash-codes to equal one another when from same Blast2GoAnnot",
				b2gaRef.hashCode() == b2gaRef.hashCode());
		assertTrue(
				"Expecting hash-code to equal Blast2GoAnnot of same description",
				b2gaRef.hashCode() == b2gaClone.hashCode());
		assertTrue(
				"Expecting hash-code NOT to equal another Blast2GoAnnot with a different description.",
				b2gaRef.hashCode() != b2gaTwo.hashCode());
	}
}
