// This code will find the highest pT jets and their particle from a ssbar and ddbar creation

#include <iostream>

#include <vector>
#include <map>
// #include "Pythia8Plugins/LHAPythia8.h" // Include the LHEF reader
#include <cstdio>
#include <sstream>  
#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include <string>
#include <fstream>
#include <algorithm> // for std::sort
#include <math.h>

using namespace Pythia8;
using namespace fastjet;
using namespace std;




double Equation_of_a_Straight_Line(double decay_distance) {
    double y_1, y_2, x_1, x_2;
    
    if (decay_distance > 5 && decay_distance <= 220) {
        y_1 = 0.95; y_2 = 0.87; x_1 = 0; x_2 = 220;
    } else if (decay_distance > 220 && decay_distance <= 380) {
        y_1 = 0.87; y_2 = 0.65; x_1 = 220; x_2 = 380;
    } else if (decay_distance > 380 && decay_distance <= 410) {
        y_1 = 0.65; y_2 = 0.2; x_1 = 380; x_2 = 410;
    } else if (decay_distance > 410 && decay_distance <= 500) {
        y_1 = 0.2; y_2 = 0.15; x_1 = 410; x_2 = 500;
    } else if (decay_distance > 0 && decay_distance <= 5) { // ONLY for charged hadrons
        y_1 = 0.95; y_2 = 0.95; x_1 = 0; x_2 = 5;
    }
    
    double a = (y_2 - y_1) / (x_2 - x_1);
    double b = y_1 - a * x_1;
    double efficiency = pow(a * decay_distance + b, 2);
    
    return efficiency;
}



int main() {

  // Creating a txt file that includes all the information you need. 
  // During the code we can add "myfile" lines and add extra information like the number of particles generated, the id or name of the final particles, etc.
  // In order to do that we should add myfile << "whatever we want to know" ;. it will add it to the text file.
  
  ofstream myfile; 
  ofstream myalphafile;

  // Set the boolean variables
  bool pT200 = false;
  bool pT45 = true;


  string filename, lheFile;

  if (pT200) {
      filename = "mymainstagger - ssbar_ddbar_creation";
      lheFile = "/Users/Edo/Simulations/MG5_aMC_v3_5_1/ssbar_ddbar_process/Events/run_01/unweighted_events.lhe";
  } else if (pT45) {
      filename = "mymainstagger - ssbar_ddbar_creation_low_energy";
      lheFile = "/Users/Edo/Simulations/MG5_aMC_v3_5_1/ssbar_ddbar_process_low_energy/Events/run_01/unweighted_events.lhe";
  }

  myfile.open(filename);
  myalphafile.open("Alpha_Development");

  Pythia pythia;
  Event& event = pythia.event;
  
  //Here we can define the number of events we want.
  //int nEvent = pythia.mode("Main:numberOfEvents = 2");
  int nEvent = 1000000;

  pythia.readString("Next:numberShowEvent = 0"); //Print event record n times
  pythia.readString("Next:numberShowInfo = 10"); //Print event information n times
  pythia.readString("Next:numberShowProcess = 10"); //Print event process record n times
  pythia.readString("Init:showChangedParticleData = 1"); //Print a list of particle and decay data for those particles that were changed (one way or another).

  pythia.readString("Beams:frameType = 4"); // Tell Pythia to read from an LHE file
  pythia.readString("Beams:LHEF = " + lheFile); // Specify the LHE file path

  // pythia.readString("Beams:idA = 2212");
  // pythia.readString("Beams:idB = 2212");
  // pythia.readString("Beams:eCM = 14000.");

  // pythia.readString("HardQCD:gg2qqbar = on");
  // pythia.readString("HardQCD:qqbar2qqbarNew  = on");
  // pythia.readString("HardQCD:qq2qq  = on");

  // pythia.readString("PhaseSpace:pThatMin = 200"); 
  // pThat is related to the momentum transfer in a collision. By setting pThatMin = 100 (for example) , we're telling Pythia that we only want to generate events where the pThat value (in units of GeV) is greater than or equal to 100. 
  //This can be useful for focusing on higher-energy events or to match the kinematic range of a certain experiment or analysis.


  // Limit decay to within a certain cylinder
  pythia.readString("ParticleDecays:limitCylinder = on");
  pythia.readString("ParticleDecays:xyMax = 1000"); // 1 meter in xy
  pythia.readString("ParticleDecays:zMax = 3000"); // 3 meter in z

  // The parameters xyMax and zMax are used to define a cylindrical volume around the origin, inside which the particles can decay.
  // This allows us to set a limit to the distance from the origin where the decay can occur. 
  // If a particle's decay vertex is located outside this cylindrical volume, the particle won't decay. 
  // This way, we can simulate the situations where the particles may decay or not decay depending on their distance from the origin.

  pythia.init();


  //______________________________________________________________________________________________________________________________________________________________//
  // Define pythia event properties:

  // Fastjet input
  // Here we creates an empty vector of PseudoJet objects that will be used to store the input particles for the jet clustering algorithm. 
  // All jets and input particles are considered PseudoJet objects.
  std::vector <fastjet::PseudoJet> fjInputs; 
  
  struct ParticleData {
    int motherIndex = 0;      // Default value of -1 (indicating not set)
    int index;            // Index of the particle in the event
    int id;               // Particle ID (310 for Ks, 3122 for Lambda)
    string name;          // Name of the particle
    double decayDistance = 0; // Decay distance from the origin in XY plane
    double productionDistance = 0; // Production distance from the origin in XY plane
    double pt;            // Transverse momentum
    double energy;        // energy
    double eta;           // Pseudorapidity
    double phi;           // Azimuthal angle
    double phi_rel;       // relative phi
    double eta_rel;       // relative eta
    bool isHadron;        // True if the particle is a hadron
    int charge;           // Charge of the particle
    bool isNeutral;       // True if the particle is neutral
    bool treatedAsNeutralHadron;
    int daughter1index = 0;
    int daughter2index = 0;
    string daughter1Name;
    string daughter2Name;
    double efficiency = 0;

  };

  
  vector<ParticleData> intermediateSParticles; // This vector holds the information (in the format defined by the struct) for each Ks and Lambda particle found in an event.
  vector<ParticleData> nonReconstructedIntermediateParticles;

  // Here we will generate the event by starting the event loop. Skip if error.
  // The next following code is into this event loop - means it is gonna run nEvent times.
  for(int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if(!pythia.next()) continue; //Generate the next event

    // Clear the intermediateSParticles and nonReconstructedIntermediateParticles vectors for the new event
    intermediateSParticles.clear();
    nonReconstructedIntermediateParticles.clear();

    // Access the process record.
    Event& process = pythia.process;

    // Access the first two partons in the process record
    // process[1] and process[2] are the protons, process[3] and process[4] are the quarks/glouns that create the process - gg/qqbar to qqbar.
    // process[5] and process[6] are the quarks we are interseted in which generated by gg/qqbar.
    Particle& parton1 = process[5]; // the first parton
    Particle& parton2 = process[6]; // the second parton

    // Check the pseudo-rapidity of the partons
    // if (fabs(parton1.eta()) > 0.05 || fabs(parton2.eta()) > 0.05) continue;

    // Get the particle ID of the first two partons. 
    int id1 = parton1.id();
    int id2 = parton2.id();

    bool isSSBar = false;
    bool isDDBar = false;

    // Check if the partons are s-sbar or d-dbar. If neither s-sbar nor d-dbar, then skip to next event.
    if ((id1 == 3 && id2 == -3) || (id1 == -3 && id2 == 3)) {
      isSSBar = true;
    } else if ((id1 == 1 && id2 == -1) || (id1 == -1 && id2 == 1)) {
      isDDBar = true;
    } else { 
      continue;
    }

    //create inputs to FastJet 
    fjInputs.resize(0);
    for (int i = 0; i < pythia.event.size(); ++i) { // Loop over all particles in the evnet 
      
      // Check if one of the intermediate particle is Ks or Lambda (Optional to add also Sigma+ 3222 and Sigma- 3112)
      if (pythia.event[i].idAbs() == 310 || pythia.event[i].idAbs() == 3122) {

          // Calculate the decay distance from the origin in the XY plane
          double decayRadiusXY = sqrt(pow(pythia.event[i].xDec(), 2) + pow(pythia.event[i].yDec(), 2));

          cout << "idx:" << i << "         "<< "name:" << pythia.event[i].name() << "         "<< "pT:" << pythia.event[i].pT() << "               "<< "Decay Radius:" << decayRadiusXY / 10 << " cm" << endl;

          // Generates a random number between 0 and 1 for each Ks particle and includes the particle only if the random number is less than 0.4.
          // Simulating a 40% efficiency in detection
          double random_number = ((double) rand() / (RAND_MAX));

          // Create a ParticleData instance and store the particle information
          ParticleData particleData;
          particleData.index = pythia.event[i].index();
          particleData.id = pythia.event[i].id();
          particleData.decayDistance = decayRadiusXY;
          particleData.pt = pythia.event[i].pT();
          particleData.energy = pythia.event[i].e();
          particleData.eta = pythia.event[i].eta();
          particleData.phi = pythia.event[i].phi();
          particleData.name = pythia.event[i].name();
          particleData.isHadron = pythia.event[i].isHadron();
          particleData.charge = pythia.event[i].charge();
          particleData.isNeutral = pythia.event[i].isNeutral();
          particleData.treatedAsNeutralHadron = false;
          
          
          particleData.daughter1index = pythia.event[i].daughter1();
          particleData.daughter2index = pythia.event[i].daughter2();
          // For daughters, first check if the daughter indices are valid
          int daughter1Index = pythia.event[i].daughter1();
          int daughter2Index = pythia.event[i].daughter2();
          
          if (daughter1Index > 0 && daughter1Index < pythia.event.size()) {
              particleData.daughter1Name = pythia.event[daughter1Index].name();
          } else {
              // particleData.daughter1Name = "No Daughter";
          }
          
          if (daughter2Index > 0 && daughter2Index < pythia.event.size()) {
              particleData.daughter2Name = pythia.event[daughter2Index].name();
          } else {
              // particleData.daughter2Name = "No Daughter";
          }


          // // Check if the decay distance is less than 50 cm and apply the efficiency condition
          // if (decayRadiusXY <= 500 && random_number > 0.4) {
          //   continue; // Skip this particle with 60% probability
          // }

          // Check if the decay distance is less than 50 cm (500 mm) and apply the efficiency condition
          if (decayRadiusXY > 5 && decayRadiusXY <= 500) {
              double efficiency = Equation_of_a_Straight_Line(decayRadiusXY);
              particleData.efficiency = efficiency;
              if (random_number <= efficiency) {

              // For a reconstructed particle (40% efficiency in detection)
              // if (random_number <= 0.4) {
                  // Add the ParticleData instance to the reconstructed vector
                  intermediateSParticles.push_back(particleData);
              }
              // For a non-reconstructed particle
              else {
                  // Add the ParticleData instance to the non-reconstructed vector
                  nonReconstructedIntermediateParticles.push_back(particleData);
              }
          } else if (decayRadiusXY > 500) {
              // For particles that decay at a distance greater than 50 cm, always add them to intermediateSParticles
              intermediateSParticles.push_back(particleData);
          }
      }

      if (!pythia.event[i].isFinal()) continue; // Accept only final stat partices
      // Explanation: if the particle is not a final state particle, the continue statement causes the loop to skip the current particle and proceed to the next particle.

      // Ignore final charged particles with pT that is lower than 0.5 and 
      if (pythia.event[i].isCharged() && pythia.event[i].pT() < 0.5) continue;
      
      // Ignore final particles with |eta| > 2.5 
      if (fabs(pythia.event[i].eta()) > 2.5) continue;


      // Ignore the neutrinos where 12, 14 and 16 are the id of the neutrinos
      if (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 || pythia.event[i].idAbs() == 16) continue;
      //cout << pythia.event[i].id() << "    " << pythia.event[i].name() << endl;  // Check the final particle that servive by printing them.
      
      // Add the particle information to the PseudoJet object 
      fastjet::PseudoJet particle = PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      particle.set_user_index(i); // Set index using set_user_info for each particle  

      fjInputs.push_back(particle); // Give fjInputs the index for each particle for later use.
      //fjInputs.push_back( fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
    }
 
    // Check that the event contains analyzable particles.
    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }
  
    // Print out the statistics
    pythia.stat();


    //______________________________________________________________________________________________________________________________________________________________//
    // Define Fastjet variabels such as R, pT miminum for each jet, the jet algorithm, etc.

    //Define the variables: anti-kt Radius, pT min
    double R = 0.4;
    // double ptmin = 5.0; //GeV 
    

    // create a jet definition: a jet algorithm with a given radius parameter
    JetDefinition jet_def(antikt_algorithm, R);

    // Run Fastjet anti-kT algorithm and sort jets in pT order.
    ClusterSequence clust_seq(fjInputs, jet_def);
    
    vector<PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets()); // cancel this line if you would like to use pT min
    // vector<PseudoJet> inclusive_jets= sorted_by_pt(clust_seq.inclusive_jets(ptmin)); // Use it if you would like to use pT min (cancel the above line)

    //Explanation:
    ///sorted_by_pt: return a vector of jets sorted into decreasing transverse momentum
    //Inclusive jets return a vector of all jets with pT > pT min
    

    //______________________________________________________________________________________________________________________________________________________________//
    // Printing the jet and constituents properties such as pT, eta and phi. In addition, we can record the number of jets per events or the constituents number per jet.


    if (inclusive_jets.size() < 2) {
      cout << "Unexpected: Less than 2 jets!" << endl;
      continue;
    }

    // print out the detils for each jet 
    for (unsigned int i = 0; i < 2 ; i++) {
    //for (unsigned int i = 0; i < inclusive_jets.size(); i++) { // Use it if you would like to use pT min
      cout << "jet number: " << i << "     "<< "pT: " << inclusive_jets[i].pt() << "     eta: " << inclusive_jets[i].eta() << "     phi: " << inclusive_jets[i].phi() << endl;
    }

    // print the jets (after we sorted them that only the jets with pT>pTmin are shown)
    for (unsigned int i = 0; i < inclusive_jets.size() ; i++) {

      if (i==2) break; // Here we choose only the 2 highest pT jets (originated from the main quarks: ssbar or ddbar)
      cout << "jet " << i << ": "<< "pT:" << inclusive_jets[i].pt() <<"   "<< "eta:" << inclusive_jets[i].eta() <<"   "<< "phi:" << inclusive_jets[i].phi() << endl; 
      
      vector<PseudoJet> constituents = inclusive_jets[i].constituents();
      
      // Here we are going to sort the constituents in a decresing pT order, from the highest to the lowest pT constituent.
      // We take the range of the constituents PseudoJet object (begining to end).
      // We define a "Lambda function" : [](const PseudoJet& c1, const PseudoJet& c2). []  - means that the lambda function has no access to any variables from the enclosing scope.
      // The lambda function takes two PseudoJet objects as input, and returns a boolean value indicating whether the first object should come before the second object in the sorted vector.
      std::sort(constituents.begin(), constituents.end(), [](const PseudoJet& c1, const PseudoJet& c2) {
        return c1.pt() > c2.pt(); // Sort by pt, in descending order
      });

      cout << "  particle #        id          name            pt             eta       eta_eff         phi       phi_eff" << endl;
  
      myfile << inclusive_jets[i].pt() << "     " << inclusive_jets[i].eta() << "     " << constituents.size();
            

      // jet_id is a variable that will get the value 0 if the jet was created form d-quark or 1 from s-quark.
      int jet_id = -1;
      
      
      








      // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
      // Save the properties of all the final constituents in 2 vectors:
      // - Neutral constituents in neutralFinalConstituents vector
      // - Other constituents in othersFinalConstituents vector

      
      vector<ParticleData> neutralFinalConstituents; // This vector will hold the neutral hadron particles.
      vector<ParticleData> othersFinalConstituents; // This vector will hold the other particles.
      
      for (unsigned j = 0; j < constituents.size(); j++) {
      
        // Get particle properties
        int index = constituents[j].user_index();

        int id = event[index].id();
        string name = event[index].name();
        bool isHadron = event[index].isHadron();
        int charge = event[index].charge();
        bool isNeutral = event[index].isNeutral();
        double pt = constituents[j].pt();
        double energy = constituents[j].e();
        double eta = constituents[j].eta();
        double phi = constituents[j].phi();

        int mother_index = event[index].mother1(); // mother1() should return the index of the mother particle

        // Calculate relative eta and phi
        double eta_rel = constituents[j].eta() - inclusive_jets[i].eta();
        double phi_rel = constituents[j].phi() - inclusive_jets[i].phi();
        if (phi_rel > M_PI) {
              phi_rel -= 2 * M_PI;
        }
        else if (phi_rel < -M_PI) {
            phi_rel += 2 * M_PI;
        }     

        // Get the production vertex coordinates
        double xProd = event[index].xProd();
        double yProd = event[index].yProd();
        // Calculate the production distance from the origin
        double productionRadiusXY = sqrt(pow(xProd, 2) + pow(yProd, 2));    


        // Create a ParticleData instance and store the particle information
        ParticleData particleData;
        particleData.motherIndex = mother_index;
        particleData.index = index;
        particleData.id = id;
        particleData.name = name;
        particleData.pt = pt;
        particleData.energy = energy;
        particleData.eta = eta;
        particleData.phi = phi;
        particleData.eta_rel = eta_rel;
        particleData.phi_rel = phi_rel;
        particleData.isHadron = isHadron;
        particleData.charge = charge;
        particleData.isNeutral = isNeutral;
        particleData.productionDistance = productionRadiusXY;

        // Add the ParticleData instance to the corresponding vector
        if (isHadron && isNeutral) {
            particleData.treatedAsNeutralHadron = true;
            neutralFinalConstituents.push_back(particleData);
        } else {
            particleData.treatedAsNeutralHadron = false;
            othersFinalConstituents.push_back(particleData);
        }
      }





      // Print intermediateSParticles before loop
      std::cout << "intermediateSParticles before loop:" << std::endl;
      for (const auto& particle : intermediateSParticles) {
          std::cout << "Index: " << particle.index << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "decayDistance: " << particle.decayDistance << "     "<< "daughter1 index: " << particle.daughter1index <<"     "<< "daughter2 index: " << particle.daughter2index << "     "<< "daughter1 name: " << particle.daughter1Name <<"     "<< "daughter2 name: " << particle.daughter2Name << std::endl;
      }

      // Print othersFinalConstituents before loop
      std::cout << "othersFinalConstituents before loop:" << std::endl;
      for (const auto& particle : othersFinalConstituents) {
          std::cout << "mother_index: " << particle.motherIndex << "     "<< "index: " << particle.index << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "productionDistance: " << particle.productionDistance << std::endl;
      }


      std::set<int> mothersAlreadyAdded; // A set to keep track of the mother particles that have already been added. This prevents duplication
      std::vector<ParticleData> particlesToAdd; // Store the new particles (mother particles) that we want to add to othersFinalConstituents
      std::vector<int> indicesToRemove; // Store the indices of the particles that we want to remove from othersFinalConstituents

      std::set<int> nonReconstructedMothersAlreadyAdded;  // A set to keep track of the non-reconstructed mother particles that have already processed. This prevents duplication

      std::cout << "Initial size of othersFinalConstituents: " << othersFinalConstituents.size() << std::endl;

      // double reconstructionEfficiency = std::sqrt(0.4);
      // // Probabilities for different cases
      // double prob_single_reconstructed = 2 * reconstructionEfficiency * (1 - reconstructionEfficiency);
      // double prob_none_reconstructed = pow((1 - reconstructionEfficiency), 2);
      // // Normalize the probabilities to sum up to 1 beacsue in the nonReconstructedIntermediateParticles vector, there are 60% of all the KS and Lambda.
      // double total = prob_single_reconstructed + prob_none_reconstructed;
      // prob_single_reconstructed = prob_single_reconstructed / total;
      // prob_none_reconstructed = prob_none_reconstructed / total;


      for (int i = 0; i < othersFinalConstituents.size(); i++) {
          int motherIndex = othersFinalConstituents[i].motherIndex; // Store the mother index of the current constituent

          // If the particle is neutral, skip to the next iteration (we don't want to look for the mother of neutron for example)
          if (othersFinalConstituents[i].isNeutral) {
              continue;
          }
          // Loop through the intermediate particles
          for (int j = 0; j < intermediateSParticles.size(); j++) {
              int intermediateIndex = intermediateSParticles[j].index;

              // If the mother index of the constituent matches the index of intermediate particle
              if (motherIndex == intermediateIndex) {

                  if (mothersAlreadyAdded.find(motherIndex) == mothersAlreadyAdded.end()) { // Check if this mother has not already been added
                      std::cout << "Adding motherIndex: " << motherIndex << " to particlesToAdd" << std::endl;
                      particlesToAdd.push_back(intermediateSParticles[j]); // Add the mother particle to particlesToAdd vector
                      mothersAlreadyAdded.insert(motherIndex); // // Mark this mother as already added
                  }
                  
                  std::cout << "Adding index " << i << " to indicesToRemove" << std::endl;
                  indicesToRemove.push_back(i);
                  break;
              }
          }

          // Loop through the non reconstructed intermediate particles
          for (int j = 0; j < nonReconstructedIntermediateParticles.size(); j++) {
              int nonReconstructedIntermediateIndex = nonReconstructedIntermediateParticles[j].index;


              double reconstructionEfficiency = sqrt(nonReconstructedIntermediateParticles[j].efficiency);
              // Probabilities for different cases
              double prob_single_reconstructed = 2 * reconstructionEfficiency * (1 - reconstructionEfficiency);
              double prob_none_reconstructed = pow((1 - reconstructionEfficiency), 2);
              // Normalize the probabilities to sum up to 1 beacsue in the random number is between 0 to 1.
              double total = prob_single_reconstructed + prob_none_reconstructed;
              prob_single_reconstructed = prob_single_reconstructed / total;
              prob_none_reconstructed = prob_none_reconstructed / total;



              // If this mother has already been processed, skip (prevents duplicating for mothers that have already been considered)
              if (nonReconstructedMothersAlreadyAdded.find(nonReconstructedIntermediateIndex) != nonReconstructedMothersAlreadyAdded.end()) {
                continue;
              }

              // If the mother index of the constituent matches the index of non reconstructed intermediate particle
              if (motherIndex == nonReconstructedIntermediateIndex) {
                // Store the indices of the daughters in the othersFinalConstituents
                std::vector<int> daughterIndices;

                // Find the two daughters in othersFinalConstituents
                for (int k = 0; k < othersFinalConstituents.size(); k++) {
                    if (othersFinalConstituents[k].motherIndex == nonReconstructedIntermediateIndex) {
                        daughterIndices.push_back(k);
                    }
                }

                // If there are two daughters, check their reconstruction
                if (daughterIndices.size() == 2) {
                    // Generate a random number to determine which daughter is reconstructed (if any)
                    double random = ((double) rand() / (RAND_MAX)); // This generates a random number between 0 and 1
                    
                    if (random <= prob_single_reconstructed) {
                        // One of the daughters is reconstructed
                        double randomDaughter = ((double) rand() / (RAND_MAX));
                        if (randomDaughter < 0.5) {
                            // First daughter is reconstructed, treat the second as a neutral hadron
                            othersFinalConstituents[daughterIndices[1]].treatedAsNeutralHadron = true;
                        } else {
                            // Second daughter is reconstructed, treat the first as a neutral hadron
                            othersFinalConstituents[daughterIndices[0]].treatedAsNeutralHadron = true;
                        }
                    } else {
                        // Neither daughter is reconstructed, treat both as neutral hadrons
                        othersFinalConstituents[daughterIndices[0]].treatedAsNeutralHadron = true;
                        othersFinalConstituents[daughterIndices[1]].treatedAsNeutralHadron = true;
                    }
                }
                nonReconstructedMothersAlreadyAdded.insert(nonReconstructedIntermediateIndex); 
                break;  
              }
          }
      }
      
      std::cout << "Number of indices to remove: " << indicesToRemove.size() << std::endl;
      // Remove indices in reverse order to avoid invalidating indices yet to be removed
      std::sort(indicesToRemove.begin(), indicesToRemove.end(), std::greater<int>());
      for (int index : indicesToRemove) {
          std::cout << "Removing index " << index << " from othersFinalConstituents" << std::endl;
          othersFinalConstituents.erase(othersFinalConstituents.begin() + index);
      } 
      
      std::cout << "Number of particles to add: " << particlesToAdd.size() << std::endl;
      // Append new particles to othersFinalConstituents
      othersFinalConstituents.insert(othersFinalConstituents.end(), particlesToAdd.begin(), particlesToAdd.end());

      std::cout << "Final size of othersFinalConstituents: " << othersFinalConstituents.size() << std::endl;
      // Print finalConstituents after sorting
      std::cout << "othersFinalConstituents before sorting:" << std::endl;
      for (const auto& particle : othersFinalConstituents) {
          std::cout << "mother_index: " << particle.motherIndex << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "productionDistance: " << particle.productionDistance << "     "<< "decayDistance: " << particle.decayDistance << std::endl;
      }


      // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
      // Transfer constituents from the othersFinalConstituents vecotor into the neutralFinalConstituents vector based on their 


      // Iterate over all the particles in the othersFinalConstituents vector
      for (auto& particle : othersFinalConstituents) {
          // If the particle is a hadron and not neutral and the production distance is greater or equal to 500
          if (particle.isHadron && !particle.isNeutral && particle.productionDistance >= 500) {
              particle.treatedAsNeutralHadron = true;
              // std::cout << "Condition 1 met by particle: "
              //     << "Index: " << particle.motherIndex << "     "
              //     << "Name: " << particle.name << "     "
              //     << "pt: " << std::setprecision(3) << particle.pt << "     "
              //     << "productionDistance: " << particle.productionDistance << "     "
              //     << "decayDistance: " << particle.decayDistance << "     "
              //     << "ID: " << particle.id << std::endl;
          }
          
                  
          // If the particle is a hadron and not neutral and the production distance is smaller to 500
          if (particle.isHadron && !particle.isNeutral && particle.productionDistance < 500) {
              double random_number_charged_hadrons = ((double) rand() / RAND_MAX);
              double efficiency_charged_hadron = Equation_of_a_Straight_Line(particle.productionDistance);
              efficiency_charged_hadron = sqrt(efficiency_charged_hadron);
              // If the random number exceeds the linear efficiency, the reconstruction fails and the particle is treated as a neutral hadron
              if (random_number_charged_hadrons > efficiency_charged_hadron) {
                 particle.treatedAsNeutralHadron = true;
              }
          }



          // If the particle is a hadron and is neutral and if the absolute id is 310 or 3122 and if the decay distance is greater or equal to 500
          if (particle.isHadron && particle.isNeutral && (abs(particle.id) == 310 || abs(particle.id) == 3122) && particle.decayDistance >= 500) {
              particle.treatedAsNeutralHadron = true;
              // std::cout << "Condition 2 met by particle: "
              //     << "Index: " << particle.motherIndex << "     "
              //     << "Name: " << particle.name << "     "
              //     << "pt: " << std::setprecision(3) << particle.pt << "     "
              //     << "productionDistance: " << particle.productionDistance << "     "
              //     << "decayDistance: " << particle.decayDistance << "     "
              //     << "ID: " << particle.id << std::endl;
          }

          // If the particle should be treated as a neutral hadron, add it to neutralFinalConstituents
          if (particle.treatedAsNeutralHadron) {
              neutralFinalConstituents.push_back(particle);
              // std::cout << "Condition 3 met by particle: "
              //     << "Index: " << particle.motherIndex << "     "
              //     << "Name: " << particle.name << "     "
              //     << "pt: " << std::setprecision(3) << particle.pt << "     "
              //     << "productionDistance: " << particle.productionDistance << "     "
              //     << "decayDistance: " << particle.decayDistance << "     "
              //     << "ID: " << particle.id << std::endl;
          }
      }




      // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
      // Building the grid for the neutral hadrons


      struct CellData {
          double ptSum = 0;
          double energySum = 0;
          double etaCenter = 0;
          double phiCenter = 0;
          int particleId = 0;
      };

      std::map<std::pair<int, int>, CellData> gridMap; //gridMpa is a map data structure, that will store the CellData of each cell using a pair of integers (eta index and phi index) as a key.


      for (auto& particle : neutralFinalConstituents) {
          
          double eta_rel = particle.eta_rel;
          double phi_rel = particle.phi_rel;
          double eta = particle.eta;
          double phi = particle.phi;
          double pt = particle.pt;
          double energy = particle.energy;

          // Calculate the indices of the cell
        //   int etaIndex = std::floor(eta_rel / 0.1);
        //   int phiIndex = std::floor(phi_rel / 0.1);

          int etaIndex = std::floor(eta / 0.1);
          int phiIndex = std::floor(phi / 0.1);

          // Get the cell from the map (it creates the cell if it doesn't exist)
          CellData &cell = gridMap[{etaIndex, phiIndex}]; // Create a new cell or check if this cell is already exist. if it, it just calls it.

          // If the particle is an sHadron one, set the property to true
        //   if ((abs(particle.id) == 310 || abs(particle.id) == 3122) && particle.decayDistance < 500) {
        //       cell.particleId = particle.id;
        //   }

          // Update the pT sum in the cell (all particles that fall into the same cell will contribute to pT Sum)
          cell.ptSum += pt;
          cell.energySum += energy;

          // Update the eta and phi center of the cell
          cell.etaCenter = (etaIndex + 0.5) * 0.1;
          cell.phiCenter = (phiIndex + 0.5) * 0.1;
          // For example, the etaIndex is 2 so it will be (2 + 0.5) * 0.1 = 0.25 which is the center of the cell from 0.2 to 0.3 in eta dimension.
      }




      // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
      // Merge the grid (neutral constituents) with the othersFinalConstituents vector (other constituents)

      vector<ParticleData> finalConstituents;

      // Insert the particles from othersFinalConstituents that have treatedAsNeutralHadron as false
      for (const auto& particle : othersFinalConstituents) {
          if (!particle.treatedAsNeutralHadron) {
              finalConstituents.push_back(particle);
          }
      }

      // Convert the cell data into ParticleData instances and add them to finalConstituents
      for (const auto& entry : gridMap) {
          const CellData& cellData = entry.second;

          // Create a ParticleData instance for this cell
          ParticleData particleData;
          particleData.pt = cellData.ptSum; // sum of pt of all neutral hadrons in this cell
          particleData.energy = cellData.energySum; // sum of energy of all neutral hadrons in this cell
          particleData.eta_rel = cellData.etaCenter; // eta center of this cell
          particleData.phi_rel = cellData.phiCenter; // phi center of this cell
          particleData.id = cellData.particleId;

          // We will set other fields to default values such as:
          particleData.index = 0; 
          particleData.name = "Cell";
          particleData.isHadron = true; 
          particleData.charge = 0;
          particleData.isNeutral = true; 

          // Add the ParticleData instance to the finalConstituents vector
          finalConstituents.push_back(particleData);
      }






      // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
      // Sort the finalConstituents vector (that includes all the constituents in the jet) by pT in descending order


      // Sort finalConstituents by pT in descending order
      std::cout << "Sorting finalConstituents" << std::endl;
      std::sort(finalConstituents.begin(), finalConstituents.end(), [](const ParticleData& p1, const ParticleData& p2) {
          return p1.pt > p2.pt; // Sort by pt, in descending order
      });
      std::cout << "Final size of finalConstituents: " << finalConstituents.size() << std::endl;

      // Print finalConstituents after sorting
      std::cout << "finalConstituents after sorting:" << std::endl;
      for (const auto& particle : finalConstituents) {
          std::cout << "mother_index: " << particle.motherIndex << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "productionDistance: " << particle.productionDistance << std::endl;
      }
      














      
      






      // Find the phi and eta of the most energetic constituent
      double phi_most_energetic = finalConstituents[0].phi_rel;
      double eta_most_energetic = finalConstituents[0].eta_rel;
      // Calculate alpha for the most energetic constituent (reference angle)
      double alpha_most_energetic = atan2(phi_most_energetic, eta_most_energetic);
    //   double alpha_most_energetic = atan(phi_most_energetic / eta_most_energetic);

    //   if (phi_most_energetic >= 0 && eta_most_energetic < 0) {
    //       alpha_most_energetic += M_PI;
    //   } else if (phi_most_energetic < 0 && eta_most_energetic < 0) {
    //       alpha_most_energetic -= M_PI;
    //   }

      // Calcualte the sum of momentum in both sides:
      double sum_momentum_right = 0; // sum of momentum fractions from alpha 0 to pi
      double sum_momentum_left = 0;  // sum of momentum fractions from alpha 0 to -pi


      // double photon_e = 0.0, electron_e = 0.0, muon_e = 0.0, hadron_e = 0.0, hadron_0_e = 0.0, ksorlambda_e = 0.0;
      double photon_e = 0.0, electron_e = 0.0, muon_e = 0.0, hadron_e = 0.0, hadron_0_e = 0.0, lambda_e = 0.0, ks_e = 0.0;
      
      // Loop through the finalConstituents vector and determine the relative pT, eta and phi in order to calcualte alpha and r.
      // for (int j = 0; j < finalConstituents.size() && j < 30; j++) {
      for (int j = 0; j < finalConstituents.size(); j++) {

          double eta_rel = finalConstituents[j].eta_rel;
          double phi_rel = finalConstituents[j].phi_rel;
   
          double pt_rel = finalConstituents[j].pt / inclusive_jets[i].pt(); // calculate the relative pT.
          // calculate eta_rel and phi_rel for each constituent with respect to the jet axis.

          // convert phi_rel and eta_rel coordinates to polar coordinates (r, alpha) where alpha is the angle relative to the most energetic constituent. 
          double r = sqrt(pow(eta_rel, 2) + pow(phi_rel, 2)); // Calculate r
          double alpha = atan2(phi_rel, eta_rel) - alpha_most_energetic; // Calculate alpha relative to the most energetic constituent
        //   double alpha = atan(phi_rel / eta_rel) - alpha_most_energetic;

        //   if (phi_rel >= 0 && eta_rel < 0) {
        //       alpha += M_PI;
        //   } else if (phi_rel < 0 && eta_rel < 0) {
        //       alpha -= M_PI;
        //   }

          //Sum the momentum fractions
          if (alpha > 0) {
               sum_momentum_right += pt_rel;
          } else if (alpha < 0) {
              sum_momentum_left += pt_rel;
          }
    

          int id = finalConstituents[j].id;
          bool isHadron = finalConstituents[j].isHadron;
          bool isNeutral = finalConstituents[j].isNeutral;
          double energy_rel = finalConstituents[j].energy / inclusive_jets[i].e(); // calculate the relative energy.
          
          
          
          // photon energy = 22
          if (id == 22) {
              photon_e += energy_rel;
          }

          // electron energy = 11
          if (id == 11 || id == -11) {
              electron_e += energy_rel;
          }

          // muon energy = 13
          if (id == 13 || id == -13) {
              muon_e += energy_rel;
          }

          // Charged hadron energy
          if (isHadron && !isNeutral) {
              hadron_e += energy_rel;
          }

          // Neutral hadron energy
          // if (isHadron && isNeutral) {
          //     if (abs(id) == 310 || abs(id) == 3122) {
          //         ksorlambda_e += energy_rel;
          //     } else if (id != 310 && id != 3122) {
          //         hadron_0_e += energy_rel;
          //     }
          // }
          if (isHadron && isNeutral) {
                if (abs(id) == 3122) {
                    lambda_e += energy_rel;
                } else if (abs(id) == 310) {
                    ks_e += energy_rel;
                } else if (id != 3122 && id != 310) {
                    hadron_0_e += energy_rel;
                }
            }
   
      }


      // myfile << "     " << photon_e << "     " << electron_e << "     " << muon_e << "     " << hadron_e << "     " << hadron_0_e << "     " << ksorlambda_e;
      myfile << "     " << photon_e << "     " << electron_e << "     " << muon_e << "     " << hadron_e << "     " << hadron_0_e << "     " << ks_e << "     " << lambda_e;



      // Decide whether to flip alpha values based on the sum of momentum fractions
      bool flip_alpha = false;
      if (sum_momentum_left > sum_momentum_right) {
          flip_alpha = true;
      }

      int photon, electron, muon, hadron, hadron_0, ksorlambda;

      // Loop through the finalConstituents vector and determine the identity of each particle, with the possibly of flipped alpha values
      // for (int j = 0; j < finalConstituents.size() && j < 30; j++) {
      for (int j = 0; j < finalConstituents.size(); j++) {

          // Get the particle ID and other properties from the finalConstituents vector
          int id = finalConstituents[j].id;
          bool isHadron = finalConstituents[j].isHadron;
          bool isNeutral = finalConstituents[j].isNeutral;
          double charge = finalConstituents[j].charge;
          double decayDistance = finalConstituents[j].decayDistance;
          double productionDistance = finalConstituents[j].productionDistance;
          double eta_rel = finalConstituents[j].eta_rel;
          double phi_rel = finalConstituents[j].phi_rel;
          double pt_rel = finalConstituents[j].pt / inclusive_jets[i].pt(); // calculate the relative pT.

          int photon = 0, electron = 0, muon = 0, hadron = 0, hadron_0 = 0, ksorlambda = 0;

          // photon = 22
          if (id == 22){
              photon = 1;
          } 

          // electron = 11
          if (id == 11){
              electron = 1;
          } else if (id == -11){
              electron = -1;
          } 

          // muon = 13
          if (id == 13){
              muon = 1;
          } else if (id == -13){
              muon = -1;
          } 


          // Charged hadron 
          if (isHadron && !isNeutral){
            if (charge > 0){
              hadron = 1;
            } else if (charge < 0){
              hadron = -1;
            } 
          } 
          
          
          // Neutral hadron
          if (isHadron && isNeutral) {
            if (abs(id) == 310 || abs(id) == 3122) {
              if (abs(id) == 310) {
                ksorlambda = -1;
              } else if (abs(id) == 3122) {
                ksorlambda = 1;
              }         
            } else if (id != 310 && id != 3122) {
                hadron_0 = 1;
            }
          } 


          // // Charged hadron 
          // if (isHadron && !isNeutral){
          //   if (productionDistance < 500){
          //     if (charge > 0){
          //       hadron = 1;
          //     } else if (charge < 0){
          //       hadron = -1;
          //     } 
          //   } else { // For productionDistance >= 500
          //     hadron_0 = 1;
          //   }
          // } 
          
          // // Neutral hadron
          // if (isHadron && isNeutral) {
          //   if (abs(id) == 310 || abs(id) == 3122) {
          //     if (decayDistance < 500) {
          //       if (abs(id) == 310) {
          //         ksorlambda = 1;
          //       } else if (abs(id) == 3122) {
          //         ksorlambda = -1;
          //       }
          //     } else {
          //         hadron_0 = 1;
          //     }
          //   } else if (id != 310 && id != 3122) {
          //       hadron_0 = 1;
          //   }
          // } 


          

          myalphafile << pt_rel << "     " << eta_rel<< "     " << phi_rel;
          
          // convert phi_rel and eta_rel coordinates to polar coordinates (r, alpha) where alpha is the angle relative to the most energetic constituent. 
          double r = sqrt(pow(eta_rel, 2) + pow(phi_rel, 2)); // Calculate r
          double alpha = atan2(phi_rel, eta_rel) - alpha_most_energetic; // Calculate alpha relative to the most energetic constituent
        //   double alpha = atan(phi_rel / eta_rel) - alpha_most_energetic;

        //   if (phi_rel >= 0 && eta_rel < 0) {
        //       alpha += M_PI;
        //   } else if (phi_rel < 0 && eta_rel < 0) {
        //       alpha -= M_PI;
        //   }

          myalphafile << "     " << alpha << "     " << r;
          
          // Flip alpha if needed
          if (flip_alpha) {
              alpha = -alpha;
          }

          myalphafile << "     " << alpha << "     ";

          myfile << "     " << photon << "     " << electron << "     " << muon << "     " << hadron << "     " << hadron_0 << "     " << ksorlambda <<  "     "  << pt_rel << "     " << r << "     " << alpha << "     ";
          
          
          cout << "   " << j <<  "               " << id<<  "         " << finalConstituents[j].name <<  "         " << finalConstituents[j].pt << "    " 
          << finalConstituents[j].eta << "    " << r << "    " << finalConstituents[j].phi << "    " << alpha << "    "  <<endl;

      }
      


      // Clear the finalConstituents list to be ready for the next jet
      finalConstituents.clear();
      mothersAlreadyAdded.clear();
      particlesToAdd.clear();
      indicesToRemove.clear();
      nonReconstructedMothersAlreadyAdded.clear();
      gridMap.clear();
    

      // Determine the label based on whether a jet was created by s-qark or d-quark
      if (isDDBar) {
        jet_id = 0;
      } else if (isSSBar) {
        jet_id = 1;
      } else {
        jet_id = 2;
      }
      myfile << jet_id << "     " << endl; //Store 0,1 or 2 whether it was a ddbar process, ssbar process or none of them found in the jet.

      myalphafile << endl;

      cout << "\n\n\n";

    }

    cout << "\n\n\n\n\n\n\n\n\n";
 
  }

  myalphafile.close();
  myfile.close();

  return 0;
}














      // struct CellData {
      //     double ptSum = 0;
      //     double etaCenter = 0;
      //     double phiCenter = 0;
      // };

      // std::map<std::pair<int, int>, CellData> gridMap; //gridMpa is a map data structure, that will store the CellData of each cell using a pair of integers (eta index and phi index) as a key.

      // vector<ParticleData> finalConstituents;      // This vector holds information for each constituent found in the jet.
      
      // for (unsigned j = 0; j < constituents.size(); j++) {
      
      //   // Get particle properties
      //   int index = constituents[j].user_index();

      //   int id = event[index].id();
      //   string name = event[index].name();
      //   bool isHadron = event[index].isHadron();
      //   int charge = event[index].charge();
      //   bool isNeutral = event[index].isNeutral();
      //   double pt = constituents[j].pt();
      //   double eta = constituents[j].eta();
      //   double phi = constituents[j].phi();

      //   int mother_index = event[index].mother1(); // mother1() should return the index of the mother particle

      //   // Calculate relative eta and phi
      //   double eta_rel = constituents[j].eta() - inclusive_jets[i].eta();
      //   double phi_rel = constituents[j].phi() - inclusive_jets[i].phi();
      //   if (phi_rel > M_PI) {
      //         phi_rel -= 2 * M_PI;
      //   }
      //   else if (phi_rel < -M_PI) {
      //       phi_rel += 2 * M_PI;
      //   }     

      //   // Get the production vertex coordinates
      //   double xProd = event[index].xProd();
      //   double yProd = event[index].yProd();
      //   // Calculate the production distance from the origin
      //   double productionRadiusXY = sqrt(pow(xProd, 2) + pow(yProd, 2));    

      //   if (isHadron && isNeutral) {

      //     // Calculate the indices of the cell
      //     int etaIndex = std::floor(eta_rel / 0.1);
      //     int phiIndex = std::floor(phi_rel / 0.1);

      //     // Get the cell from the map (it creates the cell if it doesn't exist)
      //     CellData &cell = gridMap[{etaIndex, phiIndex}]; // Create a new cell or check if this cell is already exist. if it, it just calls it.

      //     // Update the pT sum in the cell (all particles that fall into the same cell will contribute to pT Sum)
      //     cell.ptSum += pt;

      //     // Update the eta and phi center of the cell
      //     cell.etaCenter = (etaIndex + 0.5) * 0.1;
      //     cell.phiCenter = (phiIndex + 0.5) * 0.1;
      //     // For example, the etaIndex is 2 so it will be (2 + 0.5) * 0.1 = 0.25 which is the center of the cell from 0.2 to 0.3 in eta dimension.
      //   }
        
      //   if (!(isHadron && isNeutral)) {

      //     // Create a ParticleData instance and store the particle information
      //     ParticleData particleData;
      //     particleData.motherIndex = mother_index;
      //     particleData.index = index;
      //     particleData.id = id;
      //     particleData.name = name;
      //     particleData.pt = pt;
      //     particleData.eta = eta;
      //     particleData.phi = phi;
      //     particleData.eta_rel = eta_rel;
      //     particleData.phi_rel = phi_rel;
      //     particleData.isHadron = isHadron;
      //     particleData.charge = charge;
      //     particleData.isNeutral = isNeutral;
      //     particleData.productionDistance = productionRadiusXY;

      //     // Add the ParticleData instance to the finalConstituents vector
      //     finalConstituents.push_back(particleData);
      //   }
      // }


      // // Convert the cell data into ParticleData instances and add them to finalConstituents
      // for (const auto& entry : gridMap) {
      //     const CellData& cellData = entry.second;

      //     // Create a ParticleData instance for this cell
      //     ParticleData particleData;
      //     particleData.pt = cellData.ptSum; // pt is the sum of pt of all neutral hadrons in this cell
      //     particleData.eta_rel = cellData.etaCenter; // eta center of this cell
      //     particleData.phi_rel = cellData.phiCenter; // phi center of this cell

      //     // We will set other fields to default values such as:
      //     particleData.id = 0; 
      //     particleData.index = 0; 
      //     particleData.name = "Cell";
      //     particleData.isHadron = true; 
      //     particleData.charge = false;
      //     particleData.isNeutral = true; 

      //     // Add the ParticleData instance to the finalConstituents vector
      //     finalConstituents.push_back(particleData);
      // }
      
      
      // // Print intermediateSParticles before loop
      // std::cout << "intermediateSParticles before loop:" << std::endl;
      // for (const auto& particle : intermediateSParticles) {
      //     std::cout << "Index: " << particle.index << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "decayDistance: " << particle.decayDistance << std::endl;
      // }

      // // Print finalConstituents before loop
      // std::cout << "finalConstituents before loop:" << std::endl;
      // for (const auto& particle : finalConstituents) {
      //     std::cout << "mother_index: " << particle.motherIndex << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "productionDistance: " << particle.productionDistance << std::endl;
      // }


      // std::set<int> mothersAlreadyAdded; // A set to keep track of the mother particles that have already been added. This prevents duplication
      // std::vector<ParticleData> particlesToAdd; // Store the new particles (mother particles) that we want to add to finalConstituents
      // std::vector<int> indicesToRemove; // Store the indices of the particles that we want to remove from finalConstituents

      // std::set<int> nonReconstructedMothersAlreadyAdded;  // A set to keep track of the non-reconstructed mother particles that have already processed. This prevents duplication

      // std::cout << "Initial size of finalConstituents: " << finalConstituents.size() << std::endl;

      // double reconstructionEfficiency = std::sqrt(0.4);
      
      // for (int i = 0; i < finalConstituents.size(); i++) {
      //     int motherIndex = finalConstituents[i].motherIndex; // Store the mother index of the current constituent

      //     // If the particle is neutral, skip to the next iteration (we don't want to look for the mother of neutron for example)
      //     if (finalConstituents[i].isNeutral) {
      //         continue;
      //     }
      //     // Loop through the intermediate particles
      //     for (int j = 0; j < intermediateSParticles.size(); j++) {
      //         int intermediateIndex = intermediateSParticles[j].index;

      //         // If the mother index of the constituent matches the index of intermediate particle
      //         if (motherIndex == intermediateIndex) {

      //             if (mothersAlreadyAdded.find(motherIndex) == mothersAlreadyAdded.end()) { // Check if this mother has not already been added
      //                 std::cout << "Adding motherIndex: " << motherIndex << " to particlesToAdd" << std::endl;
      //                 particlesToAdd.push_back(intermediateSParticles[j]); // Add the mother particle to particlesToAdd vector
      //                 mothersAlreadyAdded.insert(motherIndex); // // Mark this mother as already added
      //             }
                  
      //             std::cout << "Adding index " << i << " to indicesToRemove" << std::endl;
      //             indicesToRemove.push_back(i);
      //             break;
      //         }
      //     }

      //     // Loop through the non reconstructed intermediate particles
      //     for (int j = 0; j < nonReconstructedIntermediateParticles.size(); j++) {
      //       int nonReconstructedIntermediateIndex = nonReconstructedIntermediateParticles[j].index;

      //       // If this mother has already been processed, skip
      //       if (nonReconstructedMothersAlreadyAdded.find(nonReconstructedIntermediateIndex) != nonReconstructedMothersAlreadyAdded.end()) {
      //         continue;
      //       }

      //       // If the mother index of the constituent matches the index of non reconstructed intermediate particle
      //       if (motherIndex == nonReconstructedIntermediateIndex) {
      //         // Store the indices of the daughters in the finalConstituents
      //         std::vector<int> daughterIndices;

      //         // Find the two daughters in finalConstituents
      //         for (int k = 0; k < finalConstituents.size(); k++) {
      //             if (finalConstituents[k].motherIndex == nonReconstructedIntermediateIndex) {
      //                 daughterIndices.push_back(k);
      //             }
      //         }

      //         // If there are two daughters, check their reconstruction
      //         if (daughterIndices.size() == 2) {
      //             // Generate a random number to determine which daughter is reconstructed (if any)
      //             double random = ((double) rand() / (RAND_MAX)); // This generates a random number between 0 and 1
                  
      //             if (random <= reconstructionEfficiency) {
      //                 // One of the daughters is reconstructed
      //                 double randomDaughter = ((double) rand() / (RAND_MAX));
      //                 if (randomDaughter < 0.5) {
      //                     // First daughter is reconstructed, treat the second as a neutral hadron
      //                     finalConstituents[daughterIndices[1]].treatedAsNeutralHadron = true;
      //                 } else {
      //                     // Second daughter is reconstructed, treat the first as a neutral hadron
      //                     finalConstituents[daughterIndices[0]].treatedAsNeutralHadron = true;
      //                 }
      //             } else {
      //                 // Neither daughter is reconstructed, treat both as neutral hadrons
      //                 finalConstituents[daughterIndices[0]].treatedAsNeutralHadron = true;
      //                 finalConstituents[daughterIndices[1]].treatedAsNeutralHadron = true;
      //             }
      //         }
      //         nonReconstructedMothersAlreadyAdded.insert(nonReconstructedIntermediateIndex); 
      //         break;  
      //       }
      //     }
      // }
      
      // std::cout << "Number of indices to remove: " << indicesToRemove.size() << std::endl;
      // // Remove indices in reverse order to avoid invalidating indices yet to be removed
      // std::sort(indicesToRemove.begin(), indicesToRemove.end(), std::greater<int>());
      // for (int index : indicesToRemove) {
      //     std::cout << "Removing index " << index << " from finalConstituents" << std::endl;
      //     finalConstituents.erase(finalConstituents.begin() + index);
      // } 
      
      // std::cout << "Number of particles to add: " << particlesToAdd.size() << std::endl;
      // // Append new particles to finalConstituents
      // finalConstituents.insert(finalConstituents.end(), particlesToAdd.begin(), particlesToAdd.end());


      // // Sort finalConstituents by pT in descending order
      // std::cout << "Sorting finalConstituents" << std::endl;
      // std::sort(finalConstituents.begin(), finalConstituents.end(), [](const ParticleData& p1, const ParticleData& p2) {
      //     return p1.pt > p2.pt; // Sort by pt, in descending order
      // });
      // std::cout << "Final size of finalConstituents: " << finalConstituents.size() << std::endl;

      // // Print finalConstituents after sorting
      // std::cout << "finalConstituents after sorting:" << std::endl;
      // for (const auto& particle : finalConstituents) {
      //     std::cout << "mother_index: " << particle.motherIndex << "     "<< "Name: " << particle.name << "     "<< "pt: " << std::setprecision(3) << particle.pt << "     "<< "productionDistance: " << particle.productionDistance << std::endl;
      // }











































// Option to iterate over 10 particles only:

     
      // if (j == 10) break;
      // myfile << "     " << id[j] << "     " <<  constituents[j].pt() << "     " << eta_eff[j] << "     " << phi_eff[j] << "     ";




// Example for a map (like gridmap):
// gridMap = {
//     {{0, 0}, {ptSum=3.2, etaCenter=0.05, phiCenter=0.05}},
//     {{0, 1}, {ptSum=1.1, etaCenter=0.05, phiCenter=0.15}},
//     {{1, 0}, {ptSum=2.8, etaCenter=0.15, phiCenter=0.05}},
// }

// With the above gridMap, the range-based for loop iterates over each entry.
// In the first iteration, entry is {{0, 0}, {ptSum=3.2, etaCenter=0.05, phiCenter=0.05}}. 
// Here, entry.first would be {0, 0}, and entry.second would be {ptSum=3.2, etaCenter=0.05, phiCenter=0.05}.
