package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.Parameters;


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
	public void testMutatePercentageBy() {
		Double mutateBy = null;
		// As it is a random variable, try a couple of times:
		for (int i = 0; i < 10000; i++) {
			mutateBy = Parameters.mutatePercentageBy();
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
		Integer mutateBy = null;
		// As it is a random variable, try a couple of times:
		for (int i = 0; i < 10000; i++) {
			mutateBy = Parameters.mutateBlastDatabaseWeightBy();
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
