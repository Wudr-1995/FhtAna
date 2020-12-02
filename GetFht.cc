#include "GetFht.h"
#include "EvtNavigator/NavBuffer.h"
#include "SniperKernel/AlgFactory.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/RecHeader.h"
#include "Event/CalibHeader.h"
#include "Event/ElecHeader.h"
#include "Event/SimHeader.h"
#include "RootWriter/RootWriter.h"

DECLARE_ALGORITHM(GetFht);

GetFht::GetFht(const std::string& name)
: AlgBase(name),
m_iEvt(0),
m_buf(0),
m_usedPmtNum(0)
{
	declProp("ChargeCut", m_qcut = 0);
	declProp("Use3inchPmt", m_3inchusedflag = false);
	declProp("Use20inchPmt", m_20inchusedflag = true);
	declProp("Pmt3inchTimeReso", m_3inchRes = 1);
	declProp("Pmt20inchTimeReso", m_20inchRes = 8);
	declProp("LSRadius", m_LSRadius = 17700);
	declProp("FileNumber", m_turn = 0);
}

bool GetFht::initialize() {
	LogDebug << "Initializing" << std::endl;
	if(not initGeomSvc())
		return false;
	SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
	if (navBuf.invalid()) {
		LogError << "Cannot get the NavBuffer @ /Event" << std::endl;
		return false;
	}
	m_buf = navBuf.data();
	// SniperPtr<RootWriter> m_rw(*getParent(), "RootWriter");
	// if (!m_rw.valid()) {
	// 	LogError << "Cannot get the RootWriter instance" << std::endl;
	// 	return false;
	// }
	// m_tree = m_rw->bookTree("Trks", "Trks");
	// m_hist = new TH2D("FhtVTheta", "FhtVTheta", )
	m_path = "/junofs/users/wudr/Muon_Sim/Real_Muon/ML/TestSample/";
    return true;
}

bool GetFht::execute() {
	LogDebug << "executing: " << m_iEvt ++ << std::endl;
	TH2D* FhtDis = new TH2D("FhtDistribution", "FhtDistribution", 100, 0, 3.14, 500, 0, 500);
	TH2D* FhtPhi = new TH2D("FhtVPhi", "FhtVPhi", 200, -3.14, 3.14, 500, 0, 500);
	// nEvt.open(m_path + "nEvt.txt");
	// int globalCount;
	// if (!nEvt.good()) {
	// 	outFile.open(m_path + "nEvt.txt");
	// 	outFile << 1 << endl;
	// 	outFile.close();
	// 	nEvt.close();
	// 	globalCount = 1;
	// }
	// else {
	// 	nEvt >> globalCount;
	// 	globalCount ++;
	// 	nEvt.close();
	// 	outFile.open(m_path + "nEvt.txt");
	// 	outFile << globalCount << endl;
	// 	outFile.close();
	// }

	JM::SimEvent* simevent = 0;
	JM::EvtNavigator* nav =m_buf->curEvt();
	std::vector<std::string>& paths = nav->getPath();
	JM::SimHeader* simheader = 0;
	for (size_t i = 0; i < paths.size(); ++i) {
		const std::string& path = paths[i];
		if (path == "/Event/SimOrig") {
			simheader = static_cast<JM::SimHeader*>(nav->getHeader("/Event/SimOrig"));
			std::cout << "SimHeader (/Event/SimOrig): " << simheader << std::endl;
			if (simheader)
				break;
		}
	}
	if (initPmt())
		LogDebug << "Initializing PMT success" << std::endl;
	else {
		LogError << "Initializing PMT fails" << std::endl;
		return true;
	}
	if (freshPmtData(FhtDis, FhtPhi))
		LogDebug << "Freshing PMT data success" << std::endl;
	else {
		LogError << "Freshing PMT data fails" << std::endl;
		delete FhtDis;
		delete FhtPhi;
		return true;
	}

	TString path = m_path + "Sample_" + m_turn + "_" + m_iEvt + ".txt";
	outFile.open(path);
	path = m_path + "SamplePhi_" + m_turn + "_" + m_iEvt + ".txt";
	LogDebug << path << std::endl;
	PhiOut.open(path);
	path = m_path + "Label_" + m_turn + "_" + m_iEvt + ".txt";
	label.open(path);
	path = m_path + "Trks_" + m_turn + "_" + m_iEvt + ".txt";
	trks.open(path);

	simevent = (JM::SimEvent*)simheader->event();
	LogDebug << "SimEventGot" << std::endl;
	nSimTrks = simevent->getTracksVec().size();
	LogDebug << "Retrieving tracks data" << std::endl;
	short NumCrossCd = 0;
	for (short i = 0; i < nSimTrks; i ++) {
		JM::SimTrack* strk = simevent->findTrackByTrkID(i);
		TVector3 Inci(strk->getInitX(), strk->getInitY(), strk->getInitZ());
		TVector3 Dir(strk->getInitPx(), strk->getInitPy(), strk->getInitPz());
		TVector3 Exit(strk->getExitX(), strk->getExitY(), strk->getExitZ());
		if (IfCrossCd(Inci, Dir, m_LSRadius)) {
			NumCrossCd ++;
			TVector3 LSInci = InciOnLS(Inci, Dir, m_LSRadius);
			TVector3 LSExit = Exit;
			if (Exit.Mag() > m_LSRadius) {
				TVector3 antiDir = - Dir;
				LSExit = InciOnLS(Exit, antiDir, m_LSRadius);
			}
			TVector3 dir = Dir.Unit();
			// trks << LSInci.X() << "\t" << LSInci.Y() << "\t" << LSInci.Z() << "\t"
			// 	<< dir.X() << "\t" << dir.Y() << "\t" << dir.Z() << endl;
			trks << LSInci.Theta() << "\t" << LSInci.Phi() << "\t"
				<< LSExit.Theta() << "\t" << LSExit.Phi() << "\t" << LSExit.Mag() << "\t"
				<< dir.Theta() << "\t" << dir.Phi() << endl;
		}
	}
	label << NumCrossCd << endl;
	for (short i = 1; i <= 100; i ++) {
		for (short j = 1; j <= 500; j ++) {
			outFile << FhtDis->GetBinContent(i, j) << "\t";
		}
		outFile << endl;
	}
	for (short i = 1; i <= 200; i ++) {
		for (short j = 1; j <= 500; j ++) {
			PhiOut << FhtPhi->GetBinContent(i, j) << "\t";
		}
		PhiOut << endl;
	}
	delete FhtDis;
	delete FhtPhi;
	outFile.close();
	PhiOut.close();
	label.close();
	trks.close();
	LogDebug << "Executed" << endl;
	return true;
}

bool GetFht::initGeomSvc() {
	SniperPtr<RecGeomSvc> rgSvc(getParent(), "RecGeomSvc");
	if (rgSvc.invalid()) {
		LogError << "Failed to get RecGeomSvc instance" << std::endl;
		return false;
	}
	m_geom = rgSvc->getCdGeom();
	return true;
}

bool GetFht::initPmt() {
	LogDebug << "Initializing PMTs" << std::endl;
	totPmtNum = 0;
	totPmtNum = m_geom->getPmtNum();
	if (!totPmtNum) {
		LogError << "Wrong PMT Number" << std::endl;
		return false;
	}
	LogDebug << "PMT Number Got" << std::endl;
	m_ptab.reserve(totPmtNum);
	m_ptab.resize(totPmtNum);
	for (unsigned int pid = 0; pid < totPmtNum; pid ++) {
		Identifier Id = Identifier(CdID::id(pid, 0));
		PmtGeom* pmt = m_geom->getPmt(Id);
		if (!pmt) {
			LogError << "Wrong PMT ID" << std::endl;
			return false;
		}
		TVector3 pmtCenter = pmt->getCenter();
		m_ptab[pid].pos = pmtCenter;
		if (CdID::is3inch(Id)) {
			m_ptab[pid].res = m_3inchRes;
			m_ptab[pid].type = _PMTINCH3;
		}
		else if (CdID::is20inch(Id)) {
			m_ptab[pid].res = m_20inchRes;
			m_ptab[pid].type = _PMTINCH20;
		}
		else {
			LogError << "Pmt[" << pid << "] is neither 3-inch or 20-inch" << std::endl;
			return false;
		}
		m_ptab[pid].q = -1;
		m_ptab[pid].fht = 99999;
		m_ptab[pid].used = false;
	}
	return true;
}

bool GetFht::freshPmtData(TH2D* ht, TH2D* pht) {
	JM::EvtNavigator* nav = m_buf->curEvt();
	if (not nav) {
		LogError << "Cannot retrieve current navigator" << std::endl;
		return false;
	}
	JM::CalibHeader* calibheader = (JM::CalibHeader*)nav->getHeader("/Event/Calib");
	if (not calibheader) {
		LogError << "Cannot retrieve '/Event/Calib'" << std::endl;
		return false;
	}
	const std::list<JM::CalibPMTChannel*>& chhlist = calibheader->event()->calibPMTCol();
	std::list<JM::CalibPMTChannel*>::const_iterator chit = chhlist.begin();
	if (chit == chhlist.end()) {
		LogDebug << "chit == chhlist.end()" << std::endl;
		return false;
	}
	while (chit != chhlist.end()) {
		JM::CalibPMTChannel* calib = *chit ++;
		Identifier id = Identifier(calib->pmtId());
		Identifier::value_type value = id.getValue();
		if (not (value & 0xFF000000) >> 24 == 0x10)
			continue;
		unsigned int pid = CdID::module(id);
		if (pid > totPmtNum) {
			LogError << "Data/Geometry Mis-Match : PmtId(" << pid << ") >= the number of PMTs." << std::endl;
			return false;
		}
		m_ptab[pid].q = calib->nPE();
		m_ptab[pid].fht = calib->firstHitTime();
		if ((CdID::is3inch(id) && m_3inchusedflag) || (CdID::is20inch(id) && m_20inchusedflag)) {
			m_ptab[pid].used = true;
			m_usedPmtNum ++;
			ht->Fill(m_ptab[pid].pos.Theta(), m_ptab[pid].fht);
			pht->Fill(m_ptab[pid].pos.Phi(), m_ptab[pid].fht);
		}
	}
	LogDebug << "Loading calibration data done" << std::endl;
	return true;
}

bool GetFht::IfCrossCd(TVector3& Inci, TVector3& Dir, Double_t R) {
	TVector3 dir = Dir.Unit();
	Double_t Dis = TMath::Sqrt(Inci.Mag() * Inci.Mag() - fabs(Inci * dir) * fabs(Inci * dir));
	LogDebug << "Distance to center: " << Dis << endl;
	if (R < Dis)
		return false;
	return true;
}

TVector3 GetFht::InciOnLS(TVector3& Inci, TVector3& Dir, Double_t R) {
	TVector3 dir = Dir.Unit();
	Double_t Dis2 = Inci.Mag() * Inci.Mag() - fabs(Inci * dir) * fabs(Inci * dir);
	TVector3 LSInci = Inci + dir * (TMath::Sqrt(Inci.Mag() * Inci.Mag() - Dis2) - TMath::Sqrt(R * R - Dis2));
	return LSInci;
}

bool GetFht::finalize() {
	LogDebug << "Finalizing" << std::endl;
	return true;
}
