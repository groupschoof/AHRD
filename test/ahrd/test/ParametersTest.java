package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.Parameters;
import ahrd.controller.Settings;

public class ParametersTest {

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testPMutateSameParameter() {
		assertEquals(0.7685474, Parameters.pMutateSameParameter(0.5), 0.000001);
		assertEquals(0, Parameters.pMutateSameParameter(0.0), 0.0);
		assertEquals(0, Parameters.pMutateSameParameter(-0.1), 0.0);
	}

	@Test
	public void testParameterToMutateRandomIndex() {
		Parameters p = getSettings().getParameters();
		Set<Integer> inds = new HashSet<Integer>();
		for (int i = 0; i < 500; i++) {
			inds.add(p.parameterToMutateRandomIndex());
		}
		assertEquals(Parameters.getNumberOfDataBaseIndependantParameters() + Parameters.getNumberOfParametersPerDataBase() * p.getBlastDatabases().size(), inds.size());
		// Parameter-Indices: 0..3 + 2 * 3 (#Blast-Databases) = 10
		// Indices 0 to 9 should all be present:
		for (int r = 0; r < 10; r++) {
			assertTrue(
					"Parameter-Index "
							+ r
							+ " should be present in Set of randomly selected Parameter-Indices.",
					inds.contains(r));
		}
	}

	@Test
	public void testNormalizeTokenScoreWeights() {
		Parameters p = getSettings().getParameters();
		p.setTokenScoreBitScoreWeight(0.5);
		p.setTokenScoreDatabaseScoreWeight(0.75);
		p.setTokenScoreOverlapScoreWeight(1.0);
		p.normalizeTokenScoreWeights();
		assertEquals(
				"The three weights in the Token-Score-Formula should sum up to 1.0",
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.0001);
	}

	@Test
	public void testMutateBlastDbWeight() {
		Parameters p = getSettings().getParameters();
		Map<String, Integer> origDbWeights = new HashMap<String, Integer>();
		for (String blastDB : p.getBlastDatabases()) {
			origDbWeights
					.put(blastDB, new Integer(p.getBlastDbWeight(blastDB)));
		}
		// test:
		for (String blastDB : p.getBlastDatabases()) {
			p.mutateBlastDatabaseWeight(blastDB);
			assertTrue(
					"Blast-Database-Weights of Blast-Database " + blastDB
							+ " should be mutated",
					!origDbWeights.get(blastDB).equals(
							p.getBlastDbWeight(blastDB)));
		}
	}

	@Test
	public void testRandomBlastDbName() {
		Parameters p = getSettings().getParameters();
		List<String> randBlastDbs = new ArrayList<String>();
		for (int i = 0; i < 1000000; i++) {
			String r = p.randomBlastDatabaseName();
			assertTrue(
					"A random Blast-Database-Name must be one of the original ones.",
					p.getBlastDatabases().contains(r));
			randBlastDbs.add(r);
		}
		for (String blastDB : p.getBlastDatabases()) {
			assertTrue(
					"Randomly selecting a Blast-Database 1,000,000 times should contain Blast-Database '"
							+ blastDB + "'.", randBlastDbs.contains(blastDB));
		}
	}

	@Test
	public void testMutateDescriptionScoreBitScoreWeight() {
		Parameters p = getSettings().getParameters();
		Map<String, Double> origDsBsWs = new HashMap<String, Double>();
		for (String blastDB : p.getBlastDatabases()) {
			origDsBsWs.put(blastDB,
					new Double(p.getDescriptionScoreBitScoreWeight(blastDB)));
		}
		// test:
		for (String blastDB : p.getBlastDatabases()) {
			p.mutateDescriptionScoreBitScoreWeight(blastDB);
			assertTrue(
					"Blast-Database-Description-Score-Bit-Score-Weights of Blast-DB "
							+ blastDB + " should have been mutated.",
					!origDsBsWs.get(blastDB).equals(
							p.getDescriptionScoreBitScoreWeight(blastDB)));
		}
	}

	@Test
	public void testMutateTokenScoreBitScoreWeight() {
		Parameters p = getSettings().getParameters();
		Double dbsw = new Double(p.getTokenScoreDatabaseScoreWeight());
		Double osw = new Double(p.getTokenScoreOverlapScoreWeight());
		Double bsw = new Double(p.getTokenScoreBitScoreWeight());
		// test:
		p.mutateTokenScoreBitScoreWeight();
		assertEquals(
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.001);
		assertTrue(
				"All three Token-Score-Weights should have changed.",
				!dbsw.equals(p.getTokenScoreDatabaseScoreWeight())
						&& !osw.equals(p.getTokenScoreOverlapScoreWeight())
						&& !bsw.equals(p.getTokenScoreBitScoreWeight()));
	}

	@Test
	public void testMutateTokenScoreOverlapScoreWeight() {
		Parameters p = getSettings().getParameters();
		Double dbsw = new Double(p.getTokenScoreDatabaseScoreWeight());
		Double osw = new Double(p.getTokenScoreOverlapScoreWeight());
		Double bsw = new Double(p.getTokenScoreBitScoreWeight());
		// test:
		p.mutateTokenScoreOverlapScoreWeight();
		assertEquals(
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.001);
		assertTrue(
				"All three Token-Score-Weights should have changed.",
				!dbsw.equals(p.getTokenScoreDatabaseScoreWeight())
						&& !osw.equals(p.getTokenScoreOverlapScoreWeight())
						&& !bsw.equals(p.getTokenScoreBitScoreWeight()));
	}

	@Test
	public void testMutateTokenScoreDatabaseScoreWeight() {
		Parameters p = getSettings().getParameters();
		Double dbsw = new Double(p.getTokenScoreDatabaseScoreWeight());
		Double osw = new Double(p.getTokenScoreOverlapScoreWeight());
		Double bsw = new Double(p.getTokenScoreBitScoreWeight());
		// test:
		p.mutateTokenScoreDatabaseScoreWeight();
		assertEquals(
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.001);
		assertTrue(
				"All three Token-Score-Weights should have changed.",
				!dbsw.equals(p.getTokenScoreDatabaseScoreWeight())
						&& !osw.equals(p.getTokenScoreOverlapScoreWeight())
						&& !bsw.equals(p.getTokenScoreBitScoreWeight()));
	}

	@Test
	public void testNeighbour() {
		// test
		Parameters n = getSettings().getParameters().neighbour(0.0);
		assertNotNull("The neighbor of current Settings mus not be NULL.", n);
		assertTrue(
				"The neighbour of current Settings must not be the same object as Settings. - Expecting a slightly changed CLONE!",
				!n.equals(getSettings().getParameters()));
		// assert "slight difference" to currently set Parameters:
		for (int i = 0; i < 500; i++) {
			n = getSettings().getParameters().neighbour(0.0);
			Settings s = getSettings();
			boolean blastParamDiff = false;
			for (String blastDbName : getSettings().getBlastDatabases()) {
				blastParamDiff = blastParamDiff
						|| (!n.getBlastDbWeight(blastDbName).equals(
								s.getBlastDbWeight(blastDbName)) || !n
								.getDescriptionScoreBitScoreWeight(blastDbName)
								.equals(s
										.getDescriptionScoreBitScoreWeight(blastDbName)));
			}
			assertTrue(
					"The cloned neighbour must differ in exactly one AHRD-parameter-field from the currently set Settings.",
					(!n.getTokenScoreBitScoreWeight().equals(s.getTokenScoreBitScoreWeight())
					|| !n.getTokenScoreDatabaseScoreWeight().equals(s.getTokenScoreDatabaseScoreWeight())
					|| !n.getTokenScoreOverlapScoreWeight().equals(s.getTokenScoreOverlapScoreWeight()))
					|| !n.getDescriptionScoreThreshold().equals(s.getDescriptionScoreThreshold())
					|| blastParamDiff);
		}	
		// Extreme Score-Increase should result in mutation of the same parameter:
		for (int paramInd = 0; paramInd < 10; paramInd++) {
			n.setLastMutatedParameter(paramInd);
			Parameters n2 = n.neighbour(1.0);
			assertEquals("The neighbor must remember which Parameter has been mutated to evolve him from its parent.",
					paramInd, n2.getLastMutatedParameter().intValue());
			if (paramInd < Parameters.getNumberOfDataBaseIndependantParameters()) {
				switch (paramInd) {
				case 0:
					assertTrue("TokenScoreBitScoreWeight should have been mutated.",
							!n.getTokenScoreBitScoreWeight().equals(n2.getTokenScoreBitScoreWeight()));
					break;
				case 1:
					assertTrue("TokenScoreDatabaseScoreWeight should have been mutated.",
							!n.getTokenScoreDatabaseScoreWeight().equals(n2.getTokenScoreDatabaseScoreWeight()));
					break;
				case 2:
					assertTrue("TokenScoreOverlapScoreWeight should have been mutated.",
							!n.getTokenScoreOverlapScoreWeight().equals(n2.getTokenScoreOverlapScoreWeight()));
					break;
				case 3:
					assertTrue("DescriptionScoreThreshold should have been mutated.",
							!n.getDescriptionScoreThreshold().equals(n2.getDescriptionScoreThreshold()));
					break;
				}
			} else {
				int indexOfBlastDbParameter = paramInd - Parameters.getNumberOfDataBaseIndependantParameters();
				String blastDbName = getSettings().getSortedBlastDatabases()
						.get((int) Math.floor((double) indexOfBlastDbParameter / Parameters.getNumberOfParametersPerDataBase()));
				switch (indexOfBlastDbParameter % Parameters.getNumberOfParametersPerDataBase()) {
				case 0:
					assertTrue("BlastDatabaseWeight of db " + blastDbName + " should have been mutated.",
							!n.getBlastDbWeight(blastDbName).equals(n2.getBlastDbWeight(blastDbName)));
					break;
				case 1:
					assertTrue("DescriptionScoreBitScoreWeight of db " + blastDbName + " should have been mutated.",
							!n.getDescriptionScoreBitScoreWeight(blastDbName).equals(n2.getDescriptionScoreBitScoreWeight(blastDbName)));
					break;
				}
			} 
		}
	}

	@Test
	public void testClone() {
		Parameters p = getSettings().getParameters();
		p.setAvgTrainingScore(0.8);
		Parameters c = p.clone();
		// test
		assertTrue("A clone should not be it's 'parent'.", p != c);
		// All primitive types get automatically cloned by Object.clone(),
		// so just test the Maps:
		for (String blastDb : p.getBlastDatabases()) {
			// For some very strange reason, calls to the getter
			// 'getBlastDbWeight(String blastDatabaseName)' return
			// Integer-Object, as in this case Integer.parseInt returns the same
			// object for the equal int-values. Maybe this is some optimization
			// of Memory-Usage in the JVM? - Anyway, this is the reason, why
			// here we compare the original Strings, held in the Map
			// blastDbSettings:
			assertTrue(
					"The BlastDatabaseWeights should not be the same objects.",
					(p.getBlastDbParameters().get(blastDb)
							.get(Settings.BLAST_DB_WEIGHT_KEY) != c
							.getBlastDbParameters().get(blastDb)
							.get(Settings.BLAST_DB_WEIGHT_KEY)));
			assertTrue(
					"The BlastDatabase's DescriptionScoreBitScoreWeight should not be the same objects.",
					System.identityHashCode(p
							.getDescriptionScoreBitScoreWeight(blastDb)) != System.identityHashCode(c
							.getDescriptionScoreBitScoreWeight(blastDb)));
		}
		// Test passing on the average training score:
		assertEquals(p.getAvgTrainingScore(), c.getAvgTrainingScore(), 0.0);
		// Assure they are different Objects. As the operator != does not reveal
		// this, we set the clone to a different value than its parent:
		c.setAvgTrainingScore(0.7);
		assertTrue(
				"Cloning should result in the average evaluation scores being different Objects.",
				p.getAvgTrainingScore() != c.getAvgTrainingScore());
	}

	@Test
	public void testEqualityAndHashCode() {
		Parameters p = getSettings().getParameters();
		// Also test evaluation-scores:
		p.setAvgTrainingScore(100.00);
		p.setAvgPrecision(0.5);
		p.setAvgRecall(0.5);
		Parameters c = p.clone();
		assertTrue("Cloned Parameters should still be equal.", p.equals(c));
		Set<Parameters> h = new HashSet<Parameters>();
		h.add(p);
		assertTrue(
				"Equal Parameters should not be added to a mathematical Set twice.",
				!h.add(c));
		// Changing scores should still fulfill equality:
		c.setAvgTrainingScore(123.45);
		assertTrue(!p.getAvgTrainingScore().equals(c.getAvgTrainingScore()));
		assertTrue("Changing scores should still fulfill equality.",
				p.equals(c));
		// Change clone and check, if it is now unequal to parent:
		c.setTokenScoreBitScoreWeight(123.00);
		assertTrue("Changed clone should not be equal to parent", !p.equals(c));
		assertTrue("Unequal Parameters should be added to a mathematical Set.",
				h.add(c));
	}

	@Test
	public void testEquals() {
		Parameters p = getSettings().getParameters();
		Parameters q = p.clone();
		assertEquals("Clones should be equal", p, q);
		q.setTokenScoreBitScoreWeight(1.5);
		assertTrue(
				"Changed TokenScoreBitScoreWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setTokenScoreDatabaseScoreWeight(1.5);
		assertTrue(
				"Changed TokenScoreDatabaseScoreWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setTokenScoreOverlapScoreWeight(1.5);
		assertTrue(
				"Changed TokenScoreOverlapScoreWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setBlastDbWeight("swissprot", "696");
		assertTrue(
				"Changed Swissprot-BlastDbWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setDescriptionScoreBitScoreWeight("swissprot", "1.5");
		assertTrue(
				"Changed Swissprot-DescriptionScoreBitScoreWeight should result in inequality.",
				!p.equals(q));
	}

	@Test
	public void testMutatePercentageBy() {
		Double mutateBy = null;
		// As it is a random variable, try a couple of times:
		for (int i = 0; i < 10000; i++) {
			mutateBy = getSettings().getParameters().mutatePercentageBy();
			// Verify:
			assertNotNull(
					"Mutation of Percentages must be done by NON NULL Double-Values.",
					mutateBy);
			assertTrue(
					"Mutation of Percentages should always be by a positive value",
					mutateBy > 0.0);
		}
	}

	@Test
	public void testMutateBlastDatabaseWeightBy() {
		getSettings().setMutatorDeviation(0.7);
		getSettings().setMutatorMean(0.7);
		Long mutateBy = null;
		// As it is a random variable, try a couple of times:
		for (int i = 0; i < 10000; i++) {
			mutateBy = getSettings().getParameters()
					.mutateBlastDatabaseWeightBy();
			// Verify:
			assertNotNull(
					"Mutation of BlastDatabaseWeights must be done by NON NULL Long-Values.",
					mutateBy);
			assertTrue(
					"Mutation of BlastDatabaseWeights should always be by a positive value, but was "
							+ mutateBy, mutateBy > 0.0);
		}
	}
}
