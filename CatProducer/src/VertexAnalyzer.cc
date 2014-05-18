#include "../interface/VertexAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	primaryVertexProducer_ = producersNames.getParameter<edm::InputTag>("primaryVertexProducer");
}

VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	primaryVertexProducer_ = producersNames.getParameter<edm::InputTag>("primaryVertexProducer");
}

VertexAnalyzer::~VertexAnalyzer()
{
}

void VertexAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootVertex)
{

	edm::Handle< reco::VertexCollection > recoVertex;
	iEvent.getByLabel(primaryVertexProducer_, recoVertex);
	if(verbosity_>1) std::cout  << "   Number of primary vertices = " << recoVertex->size() << "   Label: " << primaryVertexProducer_.label() << "   Instance: " << primaryVertexProducer_.instance() << std::endl;

	for (unsigned int j=0; j<recoVertex->size(); j++)
	{
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] Vertex  -  (vx,vy,vz)=(" << setw(12) << (*recoVertex)[j].x() << ", " << setw(12) <<  (*recoVertex)[j].y() << ", " << setw(12) << (*recoVertex)[j].z() << ")" << endl;

//		rootEvent->addPrimaryVertex( TVector3( (*recoVertex)[j].x(), (*recoVertex)[j].y(), (*recoVertex)[j].z() ) );

		const reco::Vertex * vertex = 0;
		vertex = &((*recoVertex)[j]);
		
		CatVertex localVertex(vertex->x(), vertex->y(), vertex->z());

		localVertex.setIsValid(vertex->isValid());
		localVertex.setIsFake(vertex->isFake());
		localVertex.setChi2(vertex->chi2());
		localVertex.setNdof(vertex->ndof());
		localVertex.setTracksSize(vertex->tracksSize());
		
		localVertex.setXError(vertex->xError());
		localVertex.setYError(vertex->yError());
		localVertex.setZError(vertex->zError());

		new ((*rootVertex)[j]) CatVertex (localVertex);
	}

}
