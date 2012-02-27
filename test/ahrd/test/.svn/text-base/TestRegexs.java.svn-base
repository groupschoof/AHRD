package ahrd.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;

public class TestRegexs {

	public static final String REGEX_LIST_PATH = "./test/resources/regex_list.txt";
	public static final String MATCH_LIST_PATH = "./test/resources/match_list.txt";

	@Test
	public void test() throws IOException {
		List<Pattern> patterns = compilePatterns(REGEX_LIST_PATH);
		List<String> toBeMatched = readTestMatches(MATCH_LIST_PATH);
		for (String matches : toBeMatched) {
			System.out.println("#############################");
			String log = matches;
			for (Pattern pattern : patterns) {
				Matcher matcher = pattern.matcher(matches);
				matches = matcher.replaceAll("");
				String matchingRegEx = pattern.toString();
				log += "\n" + matchingRegEx + " -> " + matches;
			}
			System.out.println(log);			
		}
	}

	private List<Pattern> compilePatterns(String pathToRegexsFile)
			throws IOException {
		List<Pattern> patternList = new ArrayList<Pattern>();
		BufferedReader br = new BufferedReader(new FileReader(new File(
				pathToRegexsFile)));
		String iterLine = null;
		while ((iterLine = br.readLine()) != null) {
			patternList.add(Pattern.compile(iterLine.trim(),
					Pattern.CASE_INSENSITIVE));
		}
		return patternList;
	}

	private List<String> readTestMatches(String pathToMatchList)
			throws IOException {
		List<String> matchesList = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(
				pathToMatchList)));
		String iterLine = null;
		while ((iterLine = br.readLine()) != null) {
			matchesList.add(iterLine.trim());
		}
		return matchesList;
	}

}
