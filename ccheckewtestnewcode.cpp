//#ifndef __CINT__

#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCEvent.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/LCTOOLS.h"
#include "EVENT/Vertex.h"
#include "IMPL/VertexImpl.h"
#include "UTIL/LCRelationNavigator.h"

//#endiF

#include <map>
#include <vector>
#include <algorithm>
#include <memory>
#include <array>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace lcio;
using namespace std;


// def my var
std::vector <int> novzerocollection; // novzero coll events
std::vector <int> vzerocollectionexist; // + vzero collection
std::vector <int> equal_pointers; // equal RecoPar pointers

std::vector <double> xmc, ymc, zmc;  // MC_positon
std::vector <double> xvzero, yvzero, zvzero;  // VZero_position
std::vector <float> diff_x, diff_y, diff_z;
std::vector <MCParticle*>  mc_vec; //  MC*_Daughters

/// -----------------------------------

#define vec_s std::vector <std::string>
#define vec std::vector <double>
#define vec_mc std::vector <std::vector <MCParticle*>>
#define vec_reco std::vector <std::vector <ReconstructedParticle*>>
#define vec_array std::vector<std::array<float,3>>  


struct my_data{	

	std::vector <std::vector <ReconstructedParticle*>> rec_vec;
	std::vector <std::vector <Vertex*>> vertex_vec;
	std::vector <std::vector <std::array<float,3>>> position;

} RECO, VZERO;


std::vector <std::array <float,3>> DIFFERENCE_ARRAY(std::array<float,3> array1, std::array<float,3> array2){
std::vector <std::array<float,3>> DIFF_VEC;	
std::array <float,3> DIFF_VAR;
		for (int j=0; j<array1.size(); j++){
				for (int i=0; i<3; i++){
					DIFF_VAR[i] = array1[i]-array2[i];}
					DIFF_VEC.push_back (DIFF_VAR);
		}

	return DIFF_VEC;
}


//#define COMPARE_VZERO struct my_data
std::vector <my_data> RECO_VEC_STRUCT;
std::vector <my_data> VZERO_VEC_STRUCT;

void ccheckewtestnewcode(const char* file) {


	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	lcReader->setReadCollectionNames( { "MCTruthRecoLink", "RecoMCTruthLink", "MCParticlesSkimmed", "V0RecoParticles", "PandoraPFOs", "V0Vertices"} );
	lcReader->open(file);

	TCanvas* c1 = new TCanvas("c1", "x MC vs RECO");
	TCanvas* c2 = new TCanvas("c2", "y MC vs RECO");
	TCanvas* c3 = new TCanvas("c3", "z MC vs RECO");
	TCanvas* c4 = new TCanvas("c4", "one");

	c1->Divide(2, 1);
	c2->Divide(2, 1);
	c3->Divide(2, 1);
	c4->Divide(2, 1);

	TH1F* histo = new TH1F("histo", "x", 100, -5, 5); // Fill 1 for no V0 Reco libraries
	TH1F* histo1 = new TH1F("histo1", "y", 100, -5, 5); // Fill 1 for no V0 Reco libraries
	TH1F* histo2 = new TH1F("histo2", "z", 100, -5, 5); // Fill 1 for no V0 Reco libraries
	TH1F* histo3 = new TH1F("histo3", "histo3", 100, -5, 5); // Fill 1 for no V0 Reco libraries
	TH1F* histo4 = new TH1F("histo4", "histo4", 100, -5, 5); // Fill 1 for no V0 Reco libraries

	TH1F* histo1vz = new TH1F("histo1vz", "x", 100, -700, 700); // Fill 1 for no V0 Reco libraries
	TH1F* histo2vz = new TH1F("histo2vz", "y", 100, -700, 700); // Fill 1 for no V0 Reco libraries
	TH1F* histo3vz = new TH1F("histo3vz", "z", 100, -700, 700); // Fill 1 for no V0 Reco libraries

	TH1F* histomc = new TH1F ("histo12", "histo12", 100, -50, 50);
	TH1F* historeco = new TH1F ("historeco", "historeco", 100, -50, 50);

	TH1F* h_diff_x  = new TH1F ("h_diff_x", "h_diff_x", 100, -300, 300);
	TH1F* h_diff_y = new TH1F ("h_diff_y", "h_diff_y", 100, -300, 300);
	TH1F* h_diff_z  = new TH1F ("h_diff_z", "h_diff_z", 100, -300, 300);

	vec_s names_histo = {"VZERO_", "MC_"};
	vec_s names_variables = {"momentum", "position"};
	vec_s coordinate = {"x_", "y_", "z_"};

	TH1F* histo_event = new TH1F ("histo_event", "histo", 100, -2000, 14800);

	int s = names_variables.size();
	for (int i = 0 ; i < s; i++ ) {
			for (auto y : coordinate) {
				y = y + names_histo.at(i) + names_variables.at(i); //xmomentum_VZERO, ymomentum_VZERo;
		}	}



	TFile* froot = new TFile( "VZeroTree.root" , "RECREATE");
	TTree* t1 = new TTree( "v0" , "reconstracted vs true");
	
	TBranch* b1 = t1->Branch("xmc", &xmc, "xmc/F");
	TBranch* b2 = t1->Branch("ymc", &ymc, "ymc/F");
	TBranch* b3 = t1->Branch("zmc", &zmc, "zmc/F");
	
	TBranch* b4 = t1->Branch("xvzero", &xvzero, "xvzero/D");
	TBranch* b5  = t1->Branch("yvzero", &yvzero, "yvzero/D");
	TBranch* b6  = t1->Branch("zvzero", &zvzero, "zvzero/D");
	

	TBranch* b7 =  t1->Branch("NO_VZERO_COLL", &novzerocollection, "novzerocollection/I");
	
	TBranch* b8 =  t1->Branch("x_diff", &diff_x,"diff_x/I");
	TBranch* b9 =  t1->Branch("y_diff", &diff_y,"diff_y/I");
   	TBranch* b10 =  t1->Branch("z_diff", &diff_z,"diff_z/I");


	int nEventsInFile = lcReader->getNumberOfEvents();
	cout << "File has " << nEventsInFile << " events" << endl;


	for (int i = 0 ; i < nEventsInFile; i++ )

	{
		
		std::cout<< "___________________ HELLO "<<std::endl;
		//-- - -  - - 
		std::vector <ReconstructedParticle*> reco_vec, asos_vec, pointer_vec;			 // RECO_MC, RECO_VZERO, RECO_VZERO_EQUAL_POINTERS
		std::vector <Vertex*> ver_reco_vec, ver_asos_vec;	
		std::vector <std::array <float,3>> pos_asos_vec;
		std::vector <std::array <float,3>>::iterator it = pos_asos_vec.begin();

		//if (i % 1000 == 0) cout << "Event " << i << endl;
		LCEvent* evt = lcReader->readNextEvent();
		LCCollection* linkCol = evt->getCollection("MCTruthRecoLink");
		LCCollection* linkRec = evt->getCollection ("RecoMCTruthLink");
		LCCollection *mcskimmed = evt->getCollection("MCParticlesSkimmed");
		LCRelationNavigator* nav = new LCRelationNavigator(linkCol);
		LCRelationNavigator* navrec = new LCRelationNavigator(linkRec);


		LCCollection *vero = nullptr;

		try {
			vero = evt -> getCollection ("V0Vertices");}
		catch (lcio::DataNotAvailableException& ) {
			cout << "its null at Vzero coll" << std::endl;}

		//MC par	
		vec_mc MC_vero;
		vec_mc daughters_vero;
		vec_mc MC_recopar;

		// VZero par
		///std::vector <Vertex*> vertex_VO
		//xvec_reco recoparticle_V0; // RECO_VEC
		//std::vector <std::array<float, 3>> position_V0; //MC_VEC

		int mn = mcskimmed->getNumberOfElements();


		for (int m = 0; m < mn; m++)
		{
			MCParticle* mc = dynamic_cast <MCParticle*> (mcskimmed->getElementAt(m));
			mc_vec.push_back(mc);
			int p = mc->getPDG();


			if  ( p == 310  || p == -310 || p == -5122 || p == 5122 || p == 22 )
			{

				std::vector <MCParticle*> daughters = mc->getDaughters();

				if (daughters.size() == 0)
				{
						continue;
									}

				double positio[3] = {mc->getVertex()[0], mc->getVertex()[1], mc->getVertex()[2]}; // rue momentum for V0 MC partlicle
				double momentu[3] = {mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2]}; // true vertex for V0 MC partlicle

				xmc.push_back(mc->getVertex()[0]); // x for daughter MC par
				ymc.push_back(mc->getVertex()[1]); // y for  daughter MC par
				zmc.push_back(mc ->getVertex()[2]); // z

				histo -> Fill(mc->getVertex()[0]);
				histo1 -> Fill(mc->getVertex()[1]);
				histo2 -> Fill(mc->getVertex()[2]);

				b1->Fill();
				b2->Fill();
				b3->Fill();

				std::cout << " ver mc is x " << positio[0] << "   y" << positio[1] << "   z" << positio [2] << "   " << std::endl;




				//_________get Weight and Recoparticles collection__________________________________________

				std::vector <float> weight;
				int max1;

				for (int j = 0; j < daughters.size(); j++)
				{
					weight = nav -> getRelatedToWeights(daughters.at(j));

					if (weight.size() == 0)
					{	continue; 	 }

					max1 = std::max_element(weight.begin(), weight.end()) - weight.begin();
					std::cout << "maxl is " << max1 << std::endl;

					for(int w = 0; w < weight.size(); w++ ) {
						ReconstructedParticle* recoparticle = dynamic_cast <ReconstructedParticle*>(nav-> getRelatedToObjects (daughters.at(j))[w]);
						//reco_vec.push_back(recoparticle);
					}



					ReconstructedParticle* recopar_mass1 = dynamic_cast <ReconstructedParticle*> (nav->getRelatedToObjects(daughters.at(j))[max1]);
					Vertex* reco_ver1 = dynamic_cast <Vertex*> (recopar_mass1->getStartVertex());
	
						if (reco_ver1 == nullptr){
						std::cout<<" null Vzero" <<std::endl;
						continue;				
					}

					std::array<float,3> pos1= { reco_ver1->getPosition()[0], reco_ver1->getPosition()[1], reco_ver1->getPosition()[2]};

					//reco_vec.push_back(recopar_mass1);
					//ver_reco_vec.push_back(reco_ver1);
					//pos_reco_vec.push_back(pos1);
					//__________________________________
				//	RECO.rec_vec = recopar_mass1;
				//	RECO.vertex_vec = reco_ver1;
					
					/*for (int i=0; i<3;i++){ 
						RECO.position[i] = pos1[i];
					};*/

				//	RECO_VEC_STRUCT.push_back(RECO);

					//____________________________________

					
					if ( vero != nullptr )

					{
					
						std::cout << std::endl;
						std::cout << "VZERO exist in Event " <<i << std::endl;
						int nVertexrec = vero->getNumberOfElements();


						for (int v = 0; v < nVertexrec; v++)
						{

							//V0 particles

							Vertex* vertexzero = dynamic_cast <Vertex*> (vero->getElementAt(v));  //  get VertexVZero Vertex for V0 (reconstracted)
							//vertex_V0.push_back(vertexzero);

							float ver_position[3] = {vertexzero->getPosition()[0], vertexzero->getPosition()[1], vertexzero->getPosition()[2]};

							histo1vz->Fill(ver_position[0]);
							histo2vz->Fill(ver_position[1]);
							histo3vz->Fill(ver_position[2]);

							xvzero.push_back( ver_position[0]);
							yvzero.push_back( ver_position[1]);
							zvzero.push_back( ver_position[2]);

							b4->Fill();
							b5->Fill();
							b6->Fill();

							//t1->Fill(); // 105228 entries

							vzerocollectionexist.push_back(i); // check in which events exist
							std::cout << " Vertex_RECO_MC is [x] " << ver_position[0] << "   [y]" << ver_position[1] << "   [z]" << ver_position [2] << "   " << std::endl;
							std::cout<<"nnnnklass"<<std::endl;

							ReconstructedParticle* asos = vertexzero->getAssociatedParticle();
							if (asos ==  nullptr){
								
								std::cout<< "NUHHHHHHHHHHHH"<<std::endl;

							}
							std::cout<< "HMMM____"<<std::endl;
							std::cout<< "type"<<asos->getType()<<std::endl;
							Vertex* ver_asos = dynamic_cast <Vertex*> (asos->getStartVertex());
							if(ver_asos==nullptr){
								std::cout<<" ITS NULL "<<std::endl;
							}

							float pos_asos[3] = {ver_asos->getPosition()[0], ver_asos->getPosition()[1],ver_asos->getPosition()[2]};
							std::cout<<"hthlgjlgdl"<<pos_asos[0]<<pos_asos[1]<<pos_asos[2]<<std::endl;

							asos_vec.push_back(asos);
							ver_asos_vec.push_back(ver_asos);
						//	std::cout<<pos_asos.size();
  						//	pos_asos_vec.insert(it, pos_asos);
							std::cout<<"LOL"<<pos_asos_vec.size()<< "size ps in loop" <<std::endl;

							//___________________________________
							
						/*	VZERO.rec_vec = asos;
							VZERO.vertex_vec = ver_asos;
							for (int b=0;b<3;b++){
								VZERO.position[b] = pos_asos[b];
							};
							VZERO_VEC_STRUCT.push_back(VZERO); */

							
							//______________________________________


							std::cout << " Type in Vertex " << asos->getType() << std::endl;
							std::cout << std::endl;

						}
						 }


					else {
					     	novzerocollection.push_back(1);
					    	histo_event-> Fill(i);
					}

				}

				//-- dot

			}

		}


	std::cout<<i<<"in event"<< "size p"<<pos_asos_vec.size()<<std::endl;

	VZERO.rec_vec.push_back(asos_vec);
	VZERO.vertex_vec.push_back(ver_asos_vec);
	VZERO.position.push_back(pos_asos_vec);

	VZERO_VEC_STRUCT.push_back(VZERO); 

}


	
// std::cout << "size of reco_vec" << std::endl;
///	std::cout << reco_vec.size() << std::endl; //122303
//	std::cout << "Size of asos" << std::endl;
//	std::cout << asos_vec.size() << std::endl; //105228

	std::vector <float> position_mc;
	std::vector <float> position_reco;
	std::vector <float> pos_diff;


/*	for ( int j = 0; j < reco_vec.size(); j++) {
		for (int l = 0; l < asos_vec.size(); l++) {

				if ( reco_vec[j] == asos_vec[l] ) 
				{
					equal_pointers.push_back(1);
					pointer_vec.push_back(reco_vec[j]);

				}
			}
		}


				
				
			////	pos_reco.push_back(pos1);
			////	pos_asos.push_back(pos2);

				float diff[3]  = { pos1[0]-pos2[0], pos1[1]-pos2[1], pos1[2]-pos2[2] };
				
				diff_x.push_back(diff[0]);
				diff_y.push_back(diff[1]);
				diff_z.push_back(diff[2]);
				
				b8->Fill();
				b9->Fill();
				b10->Fill();

			}
		}
	}
*/



//	std::cout<<"size struct"<<VZERO_VEC_STRUCT.size()<<std::endl;
	std::cout<<diff_x.size()<<std::endl; 
	std::cout << "Equal pointers - " << equal_pointers.size() << std::endl;

	c1->cd(1);
	histo->Draw();

	c1->cd(2);
	histo1vz->Draw();

	c2->cd(1);
	histo1->Draw();
	c2->cd(2);
	histo2vz->Draw();

	c3->cd(1);
	histo2->Draw();
	c3->cd(2);
	histo3vz->Draw();

	c4->cd(1);
	h_diff_x->Draw();
	c4->cd(2);
	h_diff_y->Draw();

	t1->Write();
	froot->Write();

	std::cout << "Size of One" << std::endl;
	std::cout << novzerocollection.size() << std::endl;

	//std::cout<<"Size of Equal"<< equal_coll.size()<<std::endl;
	lcReader->close();



}









