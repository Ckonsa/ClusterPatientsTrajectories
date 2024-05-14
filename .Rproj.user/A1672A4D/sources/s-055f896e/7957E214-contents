#include<Rcpp.h>
#include<map>
#include<string>
#include<limits>
using namespace Rcpp;

// Add two CharacterVectors.
// Input:
//   - v1: CharacterVector
//   - v2: CharacterVector
// Output:
//   - CharacterVector that has all the elements of the given vectors. Elements in v1 come first.
// Example:
//   - addTwoVectors(c("A", "B", "C"), c("D", "E", "F")) => c("A", "B", "C", "D", "E", "F")
//   - addTwoVectors(c("D", "E", "F"), c("A", "B", "C")) => c(D", "E", "F", "A", "B", "C")
CharacterVector addTwoVectors(CharacterVector v1, CharacterVector v2) {
  for (int i = 0; i < v2.size(); i++) {
    v1.push_back(v2[i]);
  }
  return v1;
}

// Collect data about patients.
// Finds drug trajectory for every patient.
// Input:
//   - data: dataframe
//   - IDs: numeric vector
//   - n: integer, n > 0 - length of drug trajectory
// Output:
//   - list with vectors named "IDs" and "drugs".
//     "IDs" hold information about patients IDs and
//     "drugs" contains n-length drug trajectory per patient
// [[Rcpp::export]]
List collectPatientData(DataFrame data, NumericVector IDs, int n) {
  Function getDrugTrajectory("get_patient_drug_trajectory");

  CharacterVector drugTrajectories = CharacterVector::create();
  for (int i = 0; i < IDs.size(); i++) {
    drugTrajectories = addTwoVectors(drugTrajectories, getDrugTrajectory(data, IDs[i], n));
  }
  return List::create(Named("IDs") = IDs, Named("drugs") = drugTrajectories);
}

// Add two strings alphabetically.
// Input:
//   - s1: string
//   - s2: string
// Output:
//   - string which contains given strings
// Example:
//   - addStrings("abc", "def") => "abcdef"
//   - addStrings("def", "abc") => "abcdef"
std::string addStrings(std::string s1, std::string s2) {
  std::string firstString = s1;
  std::string lastString = s2;
  if (s2 < s1) {
    firstString = s2;
    lastString = s1;
  }
  return(firstString + lastString);
}

// Compare patients drug trajectories.
// Input:
//   - patientsData: list with patient IDs and their drug trajectories
//   - punishments: numeric vector
// Output:
//   - numeric matrix with distances between patients drug trajectories
// [[Rcpp::export]]
NumericMatrix compareDrugTrajectories(List patientsData, NumericVector punishments) {
  Function distanceBetweenCodes("distance_between_codes");

  NumericVector patientsIDs = patientsData["IDs"];
  CharacterVector patientsDrugs = patientsData["drugs"];

  int matrixSize = patientsIDs.size();
  NumericMatrix matrix(matrixSize, matrixSize);

  int drugsPerPatient = patientsDrugs.size() / patientsIDs.size();
  std::map<std::string, double> drugDistances;
  for (int i = 0; i < matrixSize - 1; i++) {
    for (int j = i + 1; j < matrixSize; j++) {
      // Adds distances between all drugs of two patient
      double trajectoryDistance = 0;
      // Compares both patients 1st drugs then 2nd drugs and so on.
      for (int n = 0; n < drugsPerPatient; n++) {
        String code1 = patientsDrugs[i * drugsPerPatient + n];
        String code2 = patientsDrugs[j * drugsPerPatient + n];
        std::string codesAsOne = addStrings(code1, code2);
        // Checks if we already have the distance between these drug codes in the map
        if (drugDistances.find(codesAsOne) == drugDistances.end()) {
          drugDistances[codesAsOne] = as<double>(distanceBetweenCodes(code1, code2, punishments));
        }
        trajectoryDistance += drugDistances[codesAsOne];
      }
      matrix(i, j) = trajectoryDistance;
      matrix(j, i) = trajectoryDistance;
    }
  }
  return matrix;
}

// Compare patients drug trajectories with DTW algorithm.
// Input:
//   - patientsData: list with patient IDs and their drug trajectories
//   - punishments: numeric vector
// Output:
//   - numeric matrix with distances between patients drug trajectories
// [[Rcpp::export]]
NumericMatrix compareDrugTrajectoriesWithDTW(List patientsData, NumericVector punishments) {
  Function distanceBetweenCodes("distance_between_codes");

  NumericVector patientsIDs = patientsData["IDs"];
  CharacterVector patientsDrugs = patientsData["drugs"];
  // Matrix that holds distances between patients drug trajectories
  int matrixSize = patientsIDs.size();
  NumericMatrix matrix(matrixSize, matrixSize);

  int drugsPerPatient = patientsDrugs.size() / patientsIDs.size();
  std::map<std::string, double> drugDistances;

  // Matrix to hold different trajectory distances for DTW method
  int DTWMatrixSize = drugsPerPatient + 1;
  NumericMatrix DTWMatrix(DTWMatrixSize, DTWMatrixSize);
  for (int i = 0; i < matrixSize - 1; i++) {
    for (int j = i + 1; j < matrixSize; j++) {
      // Matrix initialization
      for (int n = 0; n < DTWMatrixSize; n++) {
        for (int m = 0; m < DTWMatrixSize; m++) {
          DTWMatrix(n, m) = std::numeric_limits<double>::infinity();
        }
      }
      DTWMatrix(0, 0) = 0;

      for (int n = 1; n < DTWMatrixSize; n++) {
        for (int m = 1; m < DTWMatrixSize; m++) {
          String code1 = patientsDrugs[i * drugsPerPatient + n - 1];
          String code2 = patientsDrugs[j * drugsPerPatient + m - 1];
          std::string codesAsOne = addStrings(code1, code2);

          double distance;
          // Checks if we already have the distance between these drug codes in the map
          if (drugDistances.find(codesAsOne) == drugDistances.end()) {
            drugDistances[codesAsOne] = as<double>(distanceBetweenCodes(code1, code2, punishments));
          }
          distance = drugDistances[codesAsOne];
          // Tries to find the minimal distance between trajectories
          DTWMatrix(n, m) = distance + std::min({DTWMatrix(n - 1, m), DTWMatrix(n, m - 1), DTWMatrix(n - 1, m - 1)});
        }
      }
      double minimalDistance = DTWMatrix(drugsPerPatient, drugsPerPatient);
      matrix(i, j) = minimalDistance;
      matrix(j, i) = minimalDistance;
    }
  }
  return matrix;
}

// Find the occurrences of given drugs per patient.
// Input:
//   - data: dataframe
//   - IDs: numeric vector
//   - drugs: character vector
// Output:
//   - numeric matrix which has count of given drugs per patient
// [[Rcpp::export]]
NumericMatrix countDrugsInCluster(DataFrame data, NumericVector IDs, CharacterVector drugs) {
  Function countAllDrugs("count_all_patient_drugs");
  Function getDrugCount("get_drug_count");

  NumericMatrix matrix(drugs.size(), IDs.size());
  for (int i = 0; i < IDs.size(); i++) {
    DataFrame drugFrequencies = countAllDrugs(data, IDs[i]);
    for (int j = 0; j < drugs.size(); j++) {
      String drug = drugs[j];
      int count = as<int>(getDrugCount(drugFrequencies, drug));
      matrix(j, i) = count;
    }
  }
  return matrix;
}

// Normalize all non-zero values in matrix with log10.
// Input:
//   - values: numeric matrix
// Output:
//   - matrix where all non-zero values have been normalized.
// [[Rcpp::export]]
NumericMatrix normalize(NumericMatrix values) {
  for (int i = 0; i < values.rows(); i++) {
    for (int j = 0; j < values.cols(); j++) {
      if (values(i, j) == 0) { continue; }
      values(i, j) = log10(values(i, j));
    }
  }
  return values;
}
