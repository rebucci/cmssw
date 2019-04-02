#include "../interface/MCExtractor.h"

// -----------------------------------------------------------------------------------------------------------
// Tokens
MCExtractor::MCExtractor(edm::EDGetTokenT< reco::GenParticleCollection > gtoken,
												 edm::EDGetTokenT< TrackingParticleCollection  > ttoken, 
                         bool doTree,
                         bool TP_hitTracker,
                         double TP_minPt,
                         double TP_maxEta,
                         double TP_maxR
                        ) 
{
  // Set everything to 0
  m_OK = false;
  
  m_gtoken = gtoken;
  m_ttoken = ttoken;
  minpt    = TP_minPt;
  maxeta   = TP_maxEta;
  maxr     = TP_maxR;

  if (TP_hitTracker) needhit = true;
  else               needhit = false;

  m_gen_x       = new std::vector<float>;
  m_gen_y       = new std::vector<float>;
  m_gen_z       = new std::vector<float>; 
  m_gen_px      = new std::vector<float>;
  m_gen_py      = new std::vector<float>;
  m_gen_pz      = new std::vector<float>;
  m_gen_proc    = new std::vector<int>;
  m_gen_pdg     = new std::vector<int>; 

  m_part_pdgId  = new std::vector<int>;
  m_part_evtId  = new std::vector<int>;
  m_part_px     = new std::vector<float>; 
  m_part_py     = new std::vector<float>; 
  m_part_pz     = new std::vector<float>;	
  m_part_pt     = new std::vector<float>;
  m_part_eta    = new std::vector<float>;   
  m_part_phi    = new std::vector<float>; 	
  m_part_x      = new std::vector<float>;   
  m_part_y      = new std::vector<float>;   
  m_part_z      = new std::vector<float>;   
  m_part_used   = new std::vector<int>;  
  m_part_stId   = new std::vector< std::vector<int> >;  
  m_hits        = new std::vector<int>;

  MCExtractor::reset();

  // Tree definition
  if (doTree) {
    m_OK = true;
    MCExtractor::createTree();
  }
}

// -----------------------------------------------------------------------------------------------------------
// Retrieval
MCExtractor::MCExtractor(TFile *a_file) {
 
  // Tree definition
  m_OK = false;

  // Variables
  m_gen_x       = new std::vector<float>;
  m_gen_y       = new std::vector<float>;
  m_gen_z       = new std::vector<float>; 
  m_gen_px      = new std::vector<float>;
  m_gen_py      = new std::vector<float>;
  m_gen_pz      = new std::vector<float>;
  m_gen_proc    = new std::vector<int>;
  m_gen_pdg     = new std::vector<int>; 

  m_part_pdgId  = new std::vector<int>;
  m_part_px     = new std::vector<float>; 
  m_part_py     = new std::vector<float>; 
  m_part_pz     = new std::vector<float>;	
  m_part_pt     = new std::vector<float>;
  m_part_eta    = new std::vector<float>;   
  m_part_phi    = new std::vector<float>; 	
  m_part_x      = new std::vector<float>;   
  m_part_y      = new std::vector<float>;   
  m_part_z      = new std::vector<float>;   
  m_part_used   = new std::vector<int>;  
  m_part_stId   = new std::vector< std::vector<int> >;  
  m_part_evtId  = new std::vector<int>;
  m_hits        = new std::vector<int>;

  MCExtractor::reset();

  // Tree
  m_tree_retrieved = dynamic_cast<TTree*>(a_file->Get("MC"));

  if (!m_tree_retrieved) {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }


  m_OK = true;

  // Branch Addresses
  m_tree_retrieved->SetBranchAddress("gen_n",   &m_gen_n);
  m_tree_retrieved->SetBranchAddress("gen_proc",&m_gen_proc);
  m_tree_retrieved->SetBranchAddress("gen_pdg", &m_gen_pdg);
  m_tree_retrieved->SetBranchAddress("gen_px",  &m_gen_px);
  m_tree_retrieved->SetBranchAddress("gen_py",  &m_gen_py);
  m_tree_retrieved->SetBranchAddress("gen_pz",  &m_gen_pz);
  m_tree_retrieved->SetBranchAddress("gen_x",   &m_gen_x);
  m_tree_retrieved->SetBranchAddress("gen_y",   &m_gen_y);
  m_tree_retrieved->SetBranchAddress("gen_z",   &m_gen_z);
  
  // Infos related to the subsequent tracking particles
  m_tree_retrieved->SetBranchAddress("subpart_n",      &m_part_n);
  m_tree_retrieved->SetBranchAddress("subpart_pdgId",  &m_part_pdgId);
  m_tree_retrieved->SetBranchAddress("subpart_evtId",  &m_part_evtId);
  m_tree_retrieved->SetBranchAddress("subpart_px",     &m_part_px);
  m_tree_retrieved->SetBranchAddress("subpart_py",     &m_part_py);
  m_tree_retrieved->SetBranchAddress("subpart_pz",     &m_part_pz);
  m_tree_retrieved->SetBranchAddress("subpart_pt",     &m_part_pt);
  m_tree_retrieved->SetBranchAddress("subpart_eta",    &m_part_eta);
  m_tree_retrieved->SetBranchAddress("subpart_phi",    &m_part_phi);
  m_tree_retrieved->SetBranchAddress("subpart_x",      &m_part_x);
  m_tree_retrieved->SetBranchAddress("subpart_y",      &m_part_y);
  m_tree_retrieved->SetBranchAddress("subpart_z",      &m_part_z);
  m_tree_retrieved->SetBranchAddress("subpart_stId",   &m_part_stId);
  m_tree_retrieved->SetBranchAddress("subpart_hits",   &m_hits);
}


// -----------------------------------------------------------------------------------------------------------
// Initialize
void MCExtractor::init(const edm::EventSetup *setup) {
  // Initializations 
  
  // Here we build the whole detector
  // We need that to retrieve all the hits
  setup->get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

}

// -----------------------------------------------------------------------------------------------------------
// Method filling the main event
void MCExtractor::writeInfo(const edm::Event *event)  { //possibly add double TP_minPt_ here?
  using namespace reco;
  
  // Reset Tree Variables :
  MCExtractor::reset();


  // First of all, we get some info on the generated event
  MCExtractor::getGenInfo(event); 


  // -----------------------------------------------------------------------------------------------------------
  // Then loop on tracking particles (TPs): 

  // Get the different Calo hits
  int n_part        = 0; // The total number of stored TPs
  
  event->getByToken(m_ttoken,TPCollection);       
  //  event->getByLabel("mix","MergedTrackTruth",TrackingParticleHandle);

  const TrackingParticleCollection tpColl = *(TPCollection.product());

  std::vector<PSimHit>::const_iterator itp; 	
      
  GlobalPoint hit_position;
  LocalPoint hitpos;

  // Loop on tracking particles 
  for (TrackingParticleCollection::size_type tpIt=0; tpIt<tpColl.size(); tpIt++) { 
    //if (n_part > m_part_nMAX) continue; // Sanity check

    TrackingParticleRef tpr(TPCollection, tpIt);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
    
    // Temporarily define tracking particle variables
    float tmp_part_pdgId = tp->pdgId();                         // Particle type
    float tmp_part_px    = tp->momentum().x();                  // >>>
    float tmp_part_py    = tp->momentum().y();                  // Momentum
    float tmp_part_pz    = tp->momentum().z();                  // <<<
    float tmp_part_eta   = tp->momentum().eta();                // Eta
    float tmp_part_phi   = tp->momentum().phi();                // Phi
    float tmp_part_x     = tp->parentVertex()->position().x();  // >>>
    float tmp_part_y     = tp->parentVertex()->position().y();  // Vertex of gen
    float tmp_part_z     = tp->parentVertex()->position().z();  // <<<
    float tmp_part_evtId = tp->eventId().rawId();  
    float tmp_hits       = tp->numberOfTrackerHits();           // Hit in Tracker?
    float tmp_part_pt    = sqrt(tmp_part_px*tmp_part_px+tmp_part_py*tmp_part_py);

    // Requirements for Tracking Particles
    if (needhit) { 
      if (tmp_hits == 0)             continue; // have hit in the tracker
    }
    if (     tmp_part_pt   < minpt)  continue; // have pT above threshold
    if (fabs(tmp_part_eta) > maxeta) continue; // have eta above threshold


    // Push Backs
    m_part_pdgId->push_back(tmp_part_pdgId);
    m_part_px   ->push_back(tmp_part_px);
    m_part_py   ->push_back(tmp_part_py);
    m_part_pz   ->push_back(tmp_part_pz);
    m_part_pt   ->push_back(tmp_part_pt);
    m_part_eta  ->push_back(tmp_part_eta);
    m_part_phi  ->push_back(tmp_part_phi);
    m_part_x    ->push_back(tmp_part_x);
    m_part_y    ->push_back(tmp_part_y);
    m_part_z    ->push_back(tmp_part_z);  
    m_part_evtId->push_back(tmp_part_evtId);  
    m_hits      ->push_back(tmp_hits);        

    the_ids.clear();

    for (TrackingParticle::g4t_iterator g4T=tp->g4Track_begin(); g4T!=tp->g4Track_end(); ++g4T) 
      the_ids.push_back(g4T->trackId());
    
    m_part_stId->push_back(the_ids); 

    ++n_part;	      

  }//end loop for on tracking particle
  
  m_part_n    = n_part;

  // Fill the tree
  MCExtractor::fillTree();
}


// -----------------------------------------------------------------------------------------------------------
// Method retrieving the generated info of the event
void MCExtractor::getGenInfo(const edm::Event *event) {
  event->getByToken(m_gtoken, genParticles);

  m_gen_n=static_cast<int>(genParticles->size());

  //std::cout << "Number of GEN particles : " << m_gen_n << std::endl; 

  int npartg=0;

  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {
    m_gen_x  ->push_back(mcIter->vx()); //
    m_gen_y  ->push_back(mcIter->vy()); // Gen of the initial MIB particle
    m_gen_z  ->push_back(mcIter->vz()); //
    m_gen_px ->push_back(mcIter->px());     //
    m_gen_py ->push_back(mcIter->py());     // Momentum
    m_gen_pz ->push_back(mcIter->pz());     //
    m_gen_pdg->push_back(mcIter->pdgId());   

    ++npartg;
  }
}


// -----------------------------------------------------------------------------------------------------------
// Method getting the info from an input file
void MCExtractor::getInfo(int ievt) {
  reset();
  m_tree_retrieved->GetEntry(ievt); 
}

// -----------------------------------------------------------------------------------------------------------
// Method initializing everything (to do before each event)
void MCExtractor::reset() {
  m_gen_n         = 0;
  m_part_n        = 0;

  m_gen_x   ->clear();
  m_gen_y   ->clear();
  m_gen_z   ->clear();
  m_gen_px  ->clear();
  m_gen_py  ->clear();  
  m_gen_pz  ->clear();  
  m_gen_proc->clear();
  m_gen_pdg ->clear(); 
  
  m_part_px   ->clear(); 
  m_part_py   ->clear(); 
  m_part_pz   ->clear(); 
  m_part_pt   ->clear();
  m_part_eta  ->clear();
  m_part_phi  ->clear();
  m_part_pdgId->clear();  
  m_part_evtId->clear();  
  m_part_stId ->clear();  
  m_part_x    ->clear();      
  m_part_y    ->clear();       
  m_part_z    ->clear();       
  m_part_used ->clear();   
  m_hits      ->clear(); 
}    

// -----------------------------------------------------------------------------------------------------------
// Fill Tree
void MCExtractor::fillTree() {
  m_tree_new->Fill(); 
}

// -----------------------------------------------------------------------------------------------------------
// Fill Size
void MCExtractor::fillSize(int size) {
  m_gen_n = size;
}

// -----------------------------------------------------------------------------------------------------------
// Get Size
int  MCExtractor::getSize() {
  return m_gen_n;
}

// -----------------------------------------------------------------------------------------------------------
// Create Tree
void MCExtractor::createTree() {
 
  m_tree_new = new TTree("MC","MC info");  

  m_tree_new->Branch("gen_n",   &m_gen_n);
  m_tree_new->Branch("gen_proc",&m_gen_proc);
  m_tree_new->Branch("gen_pdg", &m_gen_pdg);
  m_tree_new->Branch("gen_px",  &m_gen_px);
  m_tree_new->Branch("gen_py",  &m_gen_py);
  m_tree_new->Branch("gen_pz",  &m_gen_pz);
  m_tree_new->Branch("gen_x",   &m_gen_x);
  m_tree_new->Branch("gen_y",   &m_gen_y);
  m_tree_new->Branch("gen_z",   &m_gen_z);

  // Infos related to the subsequent tracking particles
  m_tree_new->Branch("subpart_n",     &m_part_n);
  m_tree_new->Branch("subpart_pdgId", &m_part_pdgId);
  m_tree_new->Branch("subpart_evtId", &m_part_evtId);
  m_tree_new->Branch("subpart_stId",  &m_part_stId);
  m_tree_new->Branch("subpart_px",    &m_part_px);
  m_tree_new->Branch("subpart_py",    &m_part_py);
  m_tree_new->Branch("subpart_pz",    &m_part_pz);
  m_tree_new->Branch("subpart_pt",    &m_part_pt);
  m_tree_new->Branch("subpart_eta",   &m_part_eta);
  m_tree_new->Branch("subpart_phi",   &m_part_phi);
  m_tree_new->Branch("subpart_x",     &m_part_x);
  m_tree_new->Branch("subpart_y",     &m_part_y);
  m_tree_new->Branch("subpart_z",     &m_part_z);
  m_tree_new->Branch("subpart_hits",  &m_hits);
} 

// -----------------------------------------------------------------------------------------------------------
// Tracking Particles
void MCExtractor::clearTP() {

  int n_TP= getNTP();

  // Loop over tracking particles
  for (int i=0;i<n_TP;++i) {
    // if (getTP_r(i)>maxr || getTP_pt(i)<minpt || fabs(getTP_eta(i))>maxeta) { //OLD
    if (getTP_r(i)>maxr || m_part_pt->at(i)<minpt || fabs(m_part_eta->at(i))>maxeta) {
      m_part_used->push_back(1);
    }
    else m_part_used->push_back(0);
  }
}

// -----------------------------------------------------------------------------------------------------------
// Matched TPs
int MCExtractor::getMatchingTP(float x,float y, float z, float px,float py, float pz) {

  int idx = -1;

  // Loop over TPs
  for (int i=0;i<getNTP();++i) { 
    if (idx!=-1) return idx;
    if (fabs(getTP_x(i)-x)>0.001)   continue;
    if (fabs(getTP_y(i)-y)>0.001)   continue;
    if (fabs(getTP_z(i)-z)>0.001)   continue;
    if (fabs(getTP_px(i)-px)>0.001) continue;
    if (fabs(getTP_py(i)-py)>0.001) continue;
    if (fabs(getTP_pz(i)-pz)>0.001) continue;
    idx=i;
  }

  return idx;
}

// -----------------------------------------------------------------------------------------------------------
// Find matched TPs
void MCExtractor::findMatchingTP(const int &stID,const int &evtID, int &itp, bool verb) {

  if (verb)
    std::cout << " Into new matching " << std::endl;

  n_TP       = getNTP();

  if (itp!=-1) return;

  // Loop over tracking particles
  for (int i=0;i<n_TP;++i) {
    if (m_part_evtId->at(i)!=evtID) continue;

    if (verb)
      std::cout << " Try to match " << std::endl;
    
    // Loop on simtrack
    for (unsigned int j=0;j<m_part_stId->at(i).size();++j) {
      if (m_part_stId->at(i).at(j)!=stID) continue;
      if (itp!=-1 && itp!=i) std::cout << " More than one tracking particle for one simtrack ID: problem!!!! " << std::endl;
      itp = i;
    }
  }

  return;
}
