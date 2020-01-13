#include "FemMembranes.hpp"
#include "Diffusion.hpp"
namespace mdx
{

  bool FemMembranes::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("FemMembranes::initialize No current mesh"));

    ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("FemMembranes::initialize No cell complex"));

    cs = &mesh->ccStructure(ccName);

    // Get the solver process
    if(!getProcess(parm("Solver Process"), solverProcess))
      throw(QString("FemMembranes::initialize Unable to make solver process: %1").arg(parm("Solver Process")));
    solverProcess->initialize(parent);
    //get the trichome process
    // if(!getProcess(parm("Trichome Process"), trichomeProcess))
    //  throw(QString("FemMembranes::initialize Unable to make trichome process: %1").arg(parm("Trichome Process")));
    //trichomeProcess->initialize(parent); 

    // Get the growth process
    if(!getProcess(parm("Growth Process"), growthProcess))
      throw(QString("FemMembranes::initialize Unable to make growth process: %1").arg(parm("Growth Process")));
    growthProcess->initialize(parent);

    // Get the subdivide process
    if(!getProcess(parm("Subdivide Process"), subdivideProcess))
      throw(QString("FemMembranes::initialize Unable to make subdivide process: %1").arg(parm("Subdivide Process")));
    subdivideProcess->initialize(parent);

    return true;
  }

  bool FemMembranes::step()
  {
    if(!solverProcess)
      throw(QString("FemMembranes::run Solver process pointer invalid"));
    if(!growthProcess)
      throw(QString("FemMembranes::run Growth process pointer invalid"));
    if(!subdivideProcess)
      throw(QString("FemMembranes::run Subdivide process pointer invalid"));

    // If converged subdivide then grow
    if(!solverProcess->step()) {

      //take a screenshot
      static int screenShotCount = 0;
      QString fileName = QString("OvuleFEMGrowth-%1.JPG").arg(screenShotCount++, 4, 10, QChar('0'));
      takeSnapshot(fileName);
 

      // assigned the trichome growth data    
      //if(!getProcess(parm("Trichome Process"), trichomeProcess))
      //  throw(QString("FemMembranes::run Unable to make trichome process: %1").arg(parm("Trichome Process")));
      //trichomeProcess->run();
      growthProcess->run();
      double growthDt = growthProcess->parm("Growth Dt").toDouble();
      growthTime += growthDt;
      if(stringToBool(solverProcess->parm("Print Stats")))
        mdxInfo << QString("Growth Step Time %1, step %2").arg(growthTime).arg(growthDt) << endl;

      // If we subdivide, we need to re-initialize the solver
      if(subdivideProcess->run()) {
        mesh->updateAll(ccName);
        // Reinitialize solver with new graph
        solverProcess->initSolver(cs);
      }
    }
    mesh->updatePositions(ccName);

    return true;
  }

  bool FemMembranes::rewind(QWidget *parent)
  {
    // To rewind, we'll reload the mesh
    Mesh *mesh = currentMesh();
    if(!mesh or mesh->file().isEmpty())
      throw(QString("No current mesh, cannot rewind"));
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    growthTime = 0;
    return meshLoad.run();
  }

  bool FemMembranes::finalize(QWidget *parent)
  {
    if(!solverProcess)
      throw(QString("FemMembranes::run Solver process pointer invalid"));

    bool result = solverProcess->finalize(parent);

    // Cleanup
    mesh = 0;
    solverProcess = 0;
    growthProcess = 0;
    subdivideProcess = 0;

    return result;
  }

  
  bool FemMembraneTrichomeProcess::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("FemMembranes::initialize No current mesh"));

    ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("FemMembranes::initialize No cell complex"));

    cs = &mesh->ccStructure(ccName);

    //double dt = femMembranes->parm("Dt");
    newSources = false;
    if(!getProcess(parm("FemMembranes Process name"), femMembranes))
      throw(QString("Can not access FemMembranes main process, check the name passed"));

    return true; 

  }
  bool FemMembraneTrichomeProcess::run()
  {
    /*Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("FemMembraneSubdivide::run No current mesh"));
    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("FemMembraneSubdivide::run Invalid cell complex"));
    //get the cell complex
    CCStructure &cs = mesh->ccStructure(ccName);
    */
    double distScale = 1000;

    //get the indexAttribute
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    //get the distance Attributes
    CCIndexDoubleAttr &hodgeAttr = mesh->attributes().attrMap<CCIndex,double>(parm("Hodge Attribute"));
    CCIndexDoubleAttr initialAttr; //= mesh->attributes().attrMap<CCIndex,double>(parm("Initial sources"));
    CCIndexDoubleAttr initialAttr2;// = mesh->attributes().attrMap<CCIndex,double>(parm("Secondary sources")); // map of the secondary sources for the secondary distance calculation (for secondary growth)
    CCIndexDoubleAttr &distanceAttr = mesh->signalAttr<double>(parm("Distance field 1 Signal"));
    CCIndexDoubleAttr &distanceAttr2 = mesh->signalAttr<double>(parm("Distance field 2 Signal"));

    //if(!getProcess("FemMembranes Process name", femMembranes))
    //  throw(QString("Can not access FemMembranes main process, check the name passed"));

    //create the two sources set for the two diffusive processes
    //clear any active selection
    MeshClearSelection(*this).run(*cs, indexAttr); 
    //load the primary Sources Set
    MeshLoadSelection(*this).run(indexAttr, parm("Initial sources Set"));  
    //initialize the sources
    InitializeFromVertexSelection(*this).run(*cs, indexAttr, initialAttr);
    
    //do the same for the secondary sources set
    MeshClearSelection(*this).run(*cs, indexAttr); 
    //load the secondary Sources Set
    MeshLoadSelection(*this).run(indexAttr, parm("Secondary sources Set"));  
    //initialize the sources
    InitializeFromVertexSelection(*this).run(*cs, indexAttr, initialAttr2);

    //check if you have sources selected at all for the distance field
    int numInitSources = 0;
    int numSecondarySources = 0;
    for(CCIndex vertex : cs->cellsOfDimension(0))
    {
       if ((initialAttr)[vertex] == 1)
         numInitSources++;

       if ((initialAttr2)[vertex] == 1)
         numSecondarySources++;

    }
    if (numInitSources == 0){
      mdxInfo << "!!!!!!!!!!!!!!!!No initial sources set for distance field is assigned, please select a source for the Distance Field to grow(under Diffusion/00 Initialize attribute from selected vertices)" << endl; 
      return true;
    }
    if (numSecondarySources == 0){
      mdxInfo << "!!!!!!!!!!!!!!!!No secondary sources set for distance field is set, please select a secondary source for the Distance Field to grow (under Diffusion/00 Initialize attribute from selected vertices) " << endl; 
      return true;
    }
    Compute2DHodgeValues(*this).run(*cs,indexAttr,hodgeAttr);
 
    // run the primary distance field calculation, for primary (main growth)
    DistanceField(*this).run(*cs,distScale,indexAttr,hodgeAttr,initialAttr,distanceAttr);
   
    //get the growth attribute  //that will be re-assigned -- this is faced-based
    fem::GrowthAttr &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute"));

    //check if the growing tricome has reached the target distance from the base
    //if that is the case, initialize 3 new sources
    //double targetDistanceFromTip = parm("Distance ring-tip").toDouble();
    double tardgeDistanceBoundaryFromTip = parm("Distance boundary-tip").toDouble();
    //double toleranceRing = parm("Tolerance ring selection").toDouble();
    double effectiveDistance = parm("Effective distance for growth").toDouble();
    double secondaryEffectiveDistance = parm("Effective distance for Secondary Growth").toDouble();

    double decayingSecondaryGrowth = parm("Decaying factor for distance Secondary Growth").toDouble();
    //double angularDistanceNodes = parm("Angular distance").toDouble();
    double scalingFactor = parm("Growth scaling factor").toDouble();
    
    //std::vector<CCIndex> ringNodes;
    std::vector<CCIndex> sourceNodes;
    double maxDistance =0;
    //double converDeg2Rad = (angularDistanceNodes * 2. * M_PI)/360.;


    //assign the primary growth field -- the main trichome growth
    for (CCIndex faces : cs->cellsOfDimension(2))
    {
      fem::Growth &fG = growthAttr[faces];
      double scaledGrowth;
      //get the average face-distance (average from vertexes)
      double faceDistance = 0;
      for(CCIndex vertex : cs->incidentCells(faces,0))
        faceDistance += (distanceAttr)[vertex];
      faceDistance *= 1./cs->bounds(faces).size();
      if (faceDistance <=  effectiveDistance) {
	double x = faceDistance / effectiveDistance;
	double x2 = x*x;
	double val = 1 + x2 * (-(22./9.) + x2 * ((17./9.) + x2 * (-(4./9.))));
	fG.kIso = scalingFactor * val;
//        fG.kIso =scalingFactor * ( 1. - (4./9. * pow((faceDistance/effectiveDistance), 6.)) + (17./9. * pow((faceDistance/effectiveDistance), 4.)) - (22./9. * pow((faceDistance/effectiveDistance), 2.) ));
      }
      else
        fG.kIso = 0.;
     }
    //look for new sources for secodnary growth, if the primary trichome has grown enough
    /// I HAVE CHANGED THIS IDEA, NODES ARE SELECTED DIRECTLY ON THE MESH //define a ring of nodes, geometrical locus of points at the same distance from the initial source
    if (!newSources){
      for(CCIndex vertex : cs->cellsOfDimension(0))
      {
        double distance = (distanceAttr)[vertex];
        if(distance > maxDistance)
          maxDistance = distance;
        //if(fabs(distance - targetDistanceFromTip) <= toleranceRing)
        //  ringNodes.push_back(vertex); 
       
      }
      mdxInfo << "maxDistance is : " << maxDistance << endl;
      //if you reach the target distance, make the nodes of the ring at the chosen angular distance sources as well
      if (maxDistance >= tardgeDistanceBoundaryFromTip)
      {  
        mdxInfo<< " I have reached maximal distance to initiate secondary growth" << endl;
        // in the next step this will tell to recalculate the distance field
        newSources = true;
      }
    }
    // if newSources are present, 
    //calculate the distance field for the secondary growth as the new sources have been defined
    // and assign the secondary growth amount 
    else
    {
      DistanceField(*this).run(*cs,distScale,indexAttr,hodgeAttr,initialAttr2,distanceAttr2);
      // assign the secondary growth field
      //first define the rescaled distance function
      //double growthTime = 0.1;
      double variableEffectiveDistance = secondaryEffectiveDistance * exp(-decayingSecondaryGrowth * femMembranes->growthTime);
      //mdxInfo << "Variable effective distance seconday growth: " << variableEffectiveDistance << endl;
      for (CCIndex faces : cs->cellsOfDimension(2))
      {
        fem::Growth &fG = growthAttr[faces];
        //get the average face-distance (average from vertexes)
        double faceDistance2 = 0;
        for(CCIndex vertex : cs->incidentCells(faces,0))
          faceDistance2 += (distanceAttr2)[vertex];
        faceDistance2 *= 1./cs->bounds(faces).size();
        //double variableEffectiveDistance = secondaryEffectiveDistance * exp(-decayingSecondaryGrowth * time)
        if (faceDistance2 <= variableEffectiveDistance ) {
	  double x = faceDistance2 / variableEffectiveDistance;
	  double x2 = x*x;
	  double val = 1 + x2 * (-(22./9.) + x2 * ((17./9.) + x2 * (-(4./9.))));
	  fG.kIso += scalingFactor * val;
        }
      }
    } 
    mesh->updatePositions(ccName);

    return true;
  };

  bool FemMembraneBisect::run()
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("FemMembraneSubdivide::run No current mesh"));
    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("FemMembraneSubdivide::run Invalid cell complex"));

    fem::ElasticTriangle3Attr &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute"));
    fem::TransIsoMaterialAttr &materialAttr = mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(parm("Material Attribute"));
    fem::PressureAttr &pressureAttr = mesh->attributes().attrMap<CCIndex, fem::Pressure>(parm("Pressure Attribute"));
    fem::PressureEdgeAttr &pressureEdgeAttr = mesh->attributes().attrMap<CCIndex, fem::PressureEdge>(parm("Pressure Edge Attribute"));
 
    fem::GrowthAttr &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute"));
    fem::DirichletAttr &dirichletAttr = mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(parm("Dirichlet Attribute"));

    fem::Triangle3Attr<Point1d> &elementMorphoAnisoAttr = mesh->attributes().attrMap<CCIndex, fem::Triangle3<Point1d>>(parm("Morphogen Anisotropy Element Attribute"));
    fem::Triangle3Attr<Point1d> &elementMorphoGrowthAttr = mesh->attributes().attrMap<CCIndex, fem::Triangle3<Point1d>>(parm("Morphogen Growth Element Attribute"));
    
    fem::MorphogenDataAttr &morphogenE2DataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Anisotropy E2 Data"));
    fem::MorphogenDataAttr &morphogenKParDataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Anisotropy KPar Data"));
    fem::MorphogenDataAttr &morphogenGrowth1DataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Growth Data 1"));
    fem::MorphogenDataAttr &morphogenGrowth2DataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Growth Data 2"));
     
    fem::MorphogenDirichletAttr &moprhogenDirichletE2Attr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Anisotropy E2 Dirichlet"));
    fem::MorphogenDirichletAttr &moprhogenDirichletKParAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Anisotropy KPar Dirichlet"));
    fem::MorphogenDirichletAttr &moprhogenDirichletGrowth1Attr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Growth Dirichlet 1"));
    fem::MorphogenDirichletAttr &moprhogenDirichletGrowth2Attr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Growth Dirichlet 2"));


    FemMembraneSubdivide sDiv(*mesh, elementAttr, elementMorphoAnisoAttr, elementMorphoGrowthAttr, materialAttr, pressureAttr, pressureEdgeAttr, dirichletAttr, growthAttr, morphogenE2DataAttr, morphogenKParDataAttr, 
                           morphogenGrowth1DataAttr, morphogenGrowth2DataAttr, moprhogenDirichletE2Attr, moprhogenDirichletKParAttr, moprhogenDirichletGrowth1Attr, moprhogenDirichletGrowth2Attr);
    mesh->updateAll(ccName);

    return run(*mesh, mesh->ccStructure(ccName), parm("Max Area").toDouble(), &sDiv);
  }



  ////////////////////////////////////////////////////////
  bool AssignCellTypeForShapeQuantifier::run()
  {
    mesh = currentMesh();
    SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &cs = mesh->ccStructure(SourceCC);
    indexAttr = &mesh->indexAttr();
    shapeAttr = &mesh->attributes().attrMap<CCIndex,CellShapeData >(parm("Shape Attribute Name"));
    
    cellType = stringToCellType(parm("Assign cell type for selected cells"));
    for(CCIndex c : cs.cellsOfDimension(3))
    {
      if ((*indexAttr)[c].selected == true){   
        (*shapeAttr)[c].cellType = cellType;
        //mdxInfo<< "Cell type assigned is: " << cellType << endl;
        //mdxInfo<< "Shape attr cell type " << (*shapeAttr)[c].cellType << endl;
      }
    }
    mesh->updateAll(SourceCC);
  };




  bool SkewSymmetricTensor::run()
  {
    mesh = currentMesh();
    SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &cs = mesh->ccStructure(SourceCC);
    indexAttr = &mesh->indexAttr();
    shapeAttr = &mesh->attributes().attrMap<CCIndex,CellShapeData >(parm("Shape Attribute Name"));
    //if(!getProcess(parm("Tissue Process"), tissueProcess))
    //  throw(QString("AuxinSolver::initialize Cannot make Tissue Process:" + parm("Tissue Process")));
    //tissueProcess->initialize(parent); 
  
    //get the barycenter, weighted on the average edge lengths for edges shared by the vtx 
    //now it is a mesh so weights are not so important...but let's use them 
    for(CCIndex c : cs.volumes()) {
      Point3d barycenter = Point3d(0., 0., 0.);
      double weightSum = 0.;
      typedef std::pair<Point3d, double> midPoint2AvArea;
      std::vector <midPoint2AvArea> midPointVec2AvArea;
      // I use the element mid point weighted by its surface for the computation
      for(CCIndex f : cs.incidentCells(c,2)){
        Point3d midPosition = Point3d(0., 0., 0.);
        double faceArea = 0;
        CCIndexVec fVertices = faceVertices(cs, f);
        midPosition = (1./3.) * ((*indexAttr)[fVertices[0]].pos + (*indexAttr)[fVertices[1]].pos + (*indexAttr)[fVertices[2]].pos);     
        faceArea = 0.5 * norm(((*indexAttr)[fVertices[1]].pos - (*indexAttr)[fVertices[0]].pos)^((*indexAttr)[fVertices[2]].pos - (*indexAttr)[fVertices[0]].pos));
        weightSum += faceArea;
        midPoint2AvArea tempPair = std::make_pair(midPosition, faceArea);
        midPointVec2AvArea.push_back(tempPair);
        barycenter += midPosition * faceArea;
      }
      //the averageLength acts as a weigth to not overextimate points which are close one to another
      barycenter *= 1./weightSum;
      (*indexAttr)[c].pos = barycenter;
      Matrix3d averagePosPos;
      double xx = 0;
      double xy = 0;
      double xz = 0;
      double yy = 0;
      double yz = 0;
      double zz = 0;
      for (uint i=0; i<midPointVec2AvArea.size(); i++){
      //for(CCIndex v : cs.incidentCells(c,0)){
        xx += midPointVec2AvArea[i].second *  pow(midPointVec2AvArea[i].first.x() - barycenter.x(), 2);

        xy += midPointVec2AvArea[i].second * (midPointVec2AvArea[i].first.x() - barycenter.x())*(midPointVec2AvArea[i].first.y() - barycenter.y());
        xz += midPointVec2AvArea[i].second * (midPointVec2AvArea[i].first.x() - barycenter.x())*(midPointVec2AvArea[i].first.z() - barycenter.z());
        yy += midPointVec2AvArea[i].second * pow(midPointVec2AvArea[i].first.y() - barycenter.y(), 2);
        yz += midPointVec2AvArea[i].second * (midPointVec2AvArea[i].first.y() - barycenter.y())*(midPointVec2AvArea[i].first.z() - barycenter.z());
        zz += midPointVec2AvArea[i].second * pow(midPointVec2AvArea[i].first.z() - barycenter.z(), 2);

      }
      xx *= 1./(weightSum);
      xy *= 1./(weightSum);
      xz *= 1./(weightSum);
      yy *= 1./(weightSum);
      yz *= 1./(weightSum);
      zz *= 1./(weightSum);
      //xx *= 1./midPointVec2AvArea.size();
      //xy *= 1./midPointVec2AvArea.size();
      //xz *= 1./midPointVec2AvArea.size();
      //yy *= 1./midPointVec2AvArea.size();
      //yz *= 1./midPointVec2AvArea.size();
      //zz *= 1./midPointVec2AvArea.size();


      averagePosPos[0] = Point3d(xx, xy, xz);
      averagePosPos[1] = Point3d(xy, yy, yz);
      averagePosPos[2] = Point3d(xz, yz, zz);
     

      Matrix3d eigVect;
      Point3d eigVal;
      eigenDecompSym3x3(averagePosPos, eigVect, eigVal);
      // the eigenvectors are already given as column vectors
      (*shapeAttr)[c].skewSymmetricTensor =  eigVect;
      (*shapeAttr)[c].skewSymmetricTensor[0][2] = eigVal[0];
      (*shapeAttr)[c].skewSymmetricTensor[1][2] = eigVal[1];
      (*shapeAttr)[c].skewSymmetricTensor[2][2] = eigVal[2];
      //(*shapeAttr)[c].skewSymmetricTensor = transpose((*shapeAttr)[c].skewSymmetricTensor);
        
      }
    return 0;
  }



  
  //compute the two anti-simmetry values wrt to the main cell anisotropy axes (as computed from the skewSymmetricTensor 
  bool AntiSymmetryTensor::run()
  {
    mesh = currentMesh();
    SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &cs = mesh->ccStructure(SourceCC);
    indexAttr = &mesh->indexAttr();
    shapeAttr = &mesh->attributes().attrMap<CCIndex,CellShapeData >(parm("Shape Attribute Name"));
    //if(!getProcess(parm("Tissue Process"), tissueProcess))
    //  throw(QString("AuxinSolver::initialize Cannot make Tissue Process:" + parm("Tissue Process")));
    //tissueProcess->initialize(parent); 
  
    //get the barycenter, weighted on the average edge lengths for edges shared by the vtx 
    //now it is a mesh so weights are not so important...but let's use them 
    for(CCIndex c : cs.cellsOfDimension(3)) {
      Point3d barycenter = Point3d(0., 0., 0.);
      double weightSum = 0.;
      typedef std::pair<Point3d, double> midPoint2AvArea;
      std::vector <midPoint2AvArea> midPointVec2AvArea;
      // I use the element mid point weighted by its surface for the computation
      for(CCIndex f : cs.incidentCells(c,2)){
        Point3d midPosition = Point3d(0., 0., 0.);
        double faceArea = 0;
        CCIndexVec fVertices = faceVertices(cs, f);
        midPosition = (1./3.) * ((*indexAttr)[fVertices[0]].pos + (*indexAttr)[fVertices[1]].pos + (*indexAttr)[fVertices[2]].pos);     
        faceArea = 0.5 * norm(((*indexAttr)[fVertices[1]].pos - (*indexAttr)[fVertices[0]].pos)^((*indexAttr)[fVertices[2]].pos - (*indexAttr)[fVertices[0]].pos));
        weightSum += faceArea;
        midPoint2AvArea tempPair = std::make_pair(midPosition, faceArea);
        midPointVec2AvArea.push_back(tempPair);
        barycenter += midPosition * faceArea;
      }
      barycenter *= 1./weightSum;
      (*indexAttr)[c].pos = barycenter;

      //rotate the coordinates so to be in the diagonal basis wrt to the skewSymmetricTensor
      
      Matrix3d rotateSkewBaisTransp;
      rotateSkewBaisTransp = transpose((*shapeAttr)[c].skewSymmetricTensor);
      Point3d minAnisoAxis = rotateSkewBaisTransp[0] ^ rotateSkewBaisTransp[1];
      rotateSkewBaisTransp[2] = minAnisoAxis;
      
      std::vector <Point3d> rotatedCellPositions;

      for(uint i=0; i< midPointVec2AvArea.size(); i++){
        rotatedCellPositions.push_back(midPointVec2AvArea[i].first);
        rotatedCellPositions[i]= rotateSkewBaisTransp* rotatedCellPositions[i];
      }
      //rotate the barycenter as well (it should not matter as the analysis is centered there..)
      //Point2d barycenter2D = Point2d(barycenter.x(), barycenter.y());
      barycenter = rotateSkewBaisTransp * barycenter;
      //get the variance  or standard deviation
      Point3d posStandDeviationSquared = Point3d(0., 0., 0.);
      for(uint i=0; i<midPointVec2AvArea.size(); i++){
         posStandDeviationSquared[0] += midPointVec2AvArea[i].second * pow((rotatedCellPositions[i].x() -barycenter.x()),2);
         posStandDeviationSquared[1] += midPointVec2AvArea[i].second * pow((rotatedCellPositions[i].y() -barycenter.y()),2);
         posStandDeviationSquared[2] += midPointVec2AvArea[i].second * pow((rotatedCellPositions[i].z() -barycenter.z()),2);

      }
      posStandDeviationSquared *= 1./weightSum;
      posStandDeviationSquared[0] = pow(posStandDeviationSquared[0], 0.5);
      posStandDeviationSquared[1] = pow(posStandDeviationSquared[1], 0.5);
      posStandDeviationSquared[2] = pow(posStandDeviationSquared[2], 0.5);

      Point3d posStandDeviation = posStandDeviationSquared;

      double asymmV1 = 0;
      double asymmV2 = 0;
      double asymmV3 = 0;
      for(uint i=0; i<midPointVec2AvArea.size(); i++){
        //get the two antisymmetry values (the first for the main anisotropy axis, the second for the min)
        asymmV1 += midPointVec2AvArea[i].second * pow((rotatedCellPositions[i].x() - barycenter.x())/posStandDeviation[0],3);
        asymmV2 += midPointVec2AvArea[i].second * pow((rotatedCellPositions[i].y() - barycenter.y())/posStandDeviation[1],3);
        asymmV3 += midPointVec2AvArea[i].second * pow((rotatedCellPositions[i].z() - barycenter.z())/posStandDeviation[2],3);

      }
      asymmV1 *= 1./weightSum;
      asymmV2 *= 1./weightSum;
      asymmV3 *= 1./weightSum;
      (*shapeAttr)[c].asymmetry = Point3d(fabs(asymmV1), fabs(asymmV2), fabs(asymmV3));
      mdxInfo << "asymm V1 is " << asymmV1 << endl;
    }
    return 0;
  }

  bool L2PericlinalSurfRatio::run()
  {
    mesh = currentMesh();
    SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &cs = mesh->ccStructure(SourceCC);
    indexAttr = &mesh->indexAttr();
    shapeAttr = &mesh->attributes().attrMap<CCIndex,CellShapeData >(parm("Shape Attribute Name"));
    spuriousContactCutoff = parm("Cutoff value for spurious contact").toDouble();
    
    mesh->updateAll(SourceCC);
    //int numL1Dome = 0;
    //int numL1 = 0;

    //create a vector of L1 and L1 Dome cells
    // and a vector of L2 cells, to save computational time
    std::vector <CCIndex> L1Cells;
    //std::vector <CCIndex> L1DomeCells;
    std::vector <CCIndex> L2Cells;

    for(CCIndex c : cs.cellsOfDimension(3)){
    //   if((*shapeAttr)[c].cellType == L1Dome)
    //     L1DomeCells.push_back(c);
       if((*shapeAttr)[c].cellType == 1)
         L1Cells.push_back(c);
         //mdxInfo<< " I am pushing back an L1cell"<< endl;
  
    }

    if (L1Cells.size() == 0 /*or L1DomeCells.size() == 0*/)
      throw(QString("60 Shape Quantifier/03 Compute L2 periclinal surface ratio::run no L1  cells selected. Plese run the process: 50 Set Cell Type and assign L1 cells"));
    
    //assign from L1 and L1 Dome neighborhood relation, L2 cells (only if they belong to the selected area).
    //also assign the L2 attribute to those cells
    //for(uint i=0; i<L1DomeCells.size(); i++){
    //  for(CCIndex c: cs.neighbors(L1DomeCells[i]))
    //  {
    //    if((*shapeAttr)[c].cellType != L1 and (*shapeAttr)[c].cellType != L1Dome)
    //      L2Cells.push_back(c);
    //  }
    //}
    std::vector <CCIndex> sporeMC;
    for(uint i=0; i<L1Cells.size(); i++){
      for(CCIndex c: cs.neighbors(L1Cells[i]))
      {
        if((*shapeAttr)[c].cellType != CellType::L1 and (*shapeAttr)[c].cellType != CellType::pSMC and (*indexAttr)[c].selected == true and (*shapeAttr)[c].cellType != CellType::CC){
          L2Cells.push_back(c);
          (*shapeAttr)[c].cellType = CellType::L2;
        }
        if((*shapeAttr)[c].cellType == CellType::pSMC)
        {
          L2Cells.push_back(c);
          sporeMC.push_back(c);
        }
        else if((*shapeAttr)[c].cellType == CellType::CC)
          L2Cells.push_back(c);
        
      }
    }
    //now we should have the full list of L2 cells, even redundant as also some I do not have in the active selection will be there.
    //for each L2 cell, if selected, find the portion of cell wall shared with L1 Dome cell and with cell not belonging to L1, L1 Dome and L2
   
    for(uint i=0; i<L2Cells.size(); i++){
     double sharedTopWallArea = 0;
     double sharedBottomWallArea = 0;
     for(CCIndex nL2: cs.neighbors(L2Cells[i]))
     {
         //if(cs.adjacent(L2Cells[i],L1Cells[j]))
         //{
         if((*shapeAttr)[nL2].cellType == L1){
           for(CCIndex f: cs.incidentCells(L2Cells[i],2))
           {
             if(cs.incident(f,nL2))
               sharedTopWallArea += (*indexAttr)[f].measure;
           }
         }
         else if((*shapeAttr)[nL2].cellType != CellType::L2 and (*shapeAttr)[nL2].cellType != CellType::pSMC and (*shapeAttr)[nL2].cellType != CellType::CC ){
           if ((*indexAttr)[nL2].selected == true){
             (*shapeAttr)[nL2].cellType = L3;
             for(CCIndex f: cs.incidentCells(L2Cells[i],2))
             {
               if(cs.incident(f,nL2))
               sharedBottomWallArea += (*indexAttr)[f].measure;
             }
           }
         }

                    
     }
     (*shapeAttr)[L2Cells[i]].bottomPericlinalWallArea = sharedBottomWallArea;
     (*shapeAttr)[L2Cells[i]].topPericlinalWallArea = sharedTopWallArea;

     if (sharedBottomWallArea > 0.){
       double ratio = sharedTopWallArea/sharedBottomWallArea;
       if (sharedBottomWallArea < spuriousContactCutoff)
         (*shapeAttr)[L2Cells[i]].L2periclinalRatio = -1;
       else
         (*shapeAttr)[L2Cells[i]].L2periclinalRatio = ratio;
     }
     else 
       (*shapeAttr)[L2Cells[i]].L2periclinalRatio = -1;
   }  
   //assign CC cells authomatically if not done by user already  
   for(uint i=0; i<sporeMC.size(); i++){
     for(CCIndex nSMC: cs.neighbors(sporeMC[i]))
     {
       if ((*shapeAttr)[nSMC].cellType == CellType::L2)
         (*shapeAttr)[nSMC].cellType = CellType::CC;
     }
   }
  };


  bool VisualizeShapeQuantifiers::run(Mesh *mesh)
  {
    //mesh = currentMesh();
    SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    OutputCC = parm("Output CC");
    AnisotropyVecSize = parm("Anisotropy Vector Size").toDouble();
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &src = mesh->ccStructure(SourceCC);
    CCStructure &out = mesh->ccStructure(OutputCC);
    indexAttr = &mesh->indexAttr();
    shapeAttr = &mesh->attributes().attrMap<CCIndex,CellShapeData >(parm("Shape Attribute Name"));
    CCIndexDoubleAttr &anisoAttr = mesh->signalAttr<double>("CellAnisotropySignal");
    CCIndexDoubleAttr &maxElongAttr = mesh->signalAttr<double>("CellElongationSignal");
    CCIndexDoubleAttr &maxAsymmetryAttr = mesh->signalAttr<double>("CellMaxAsymmetrySignal");
    CCIndexDoubleAttr &midAsymmetryAttr = mesh->signalAttr<double>("CellMidAsymmetrySignal");
    CCIndexDoubleAttr &minAsymmetryAttr = mesh->signalAttr<double>("CellMinAsymmetrySignal");  
    CCIndexDoubleAttr &L2PericlinalRatioAttr = mesh->signalAttr<double>("L2PericlinalRation"); 
    CCIndexDoubleAttr &CellTypeAttr = mesh->signalAttr<double>("CellTypeAttr");      
    out = CCStructure(2);
    forall(CCIndex f, src.cellsOfDimension(3)) {
       //as a scalar I plot the max Dimension/ mid Dimension
       (anisoAttr)[f] = ((*shapeAttr)[f].skewSymmetricTensor[0][2] / (*shapeAttr)[f].skewSymmetricTensor[1][2]);  
       (maxElongAttr)[f]= (*shapeAttr)[f].skewSymmetricTensor[0][2];
       (maxAsymmetryAttr)[f] = (*shapeAttr)[f].asymmetry[0];
       (midAsymmetryAttr)[f] = (*shapeAttr)[f].asymmetry[1];
       (minAsymmetryAttr)[f] = (*shapeAttr)[f].asymmetry[2];
       (L2PericlinalRatioAttr)[f] = (*shapeAttr)[f].L2periclinalRatio;
       (CellTypeAttr)[f] = (*shapeAttr)[f].cellType;
       //if(DrawPolarizer) {
       CCIndexData V = (*indexAttr)[f];
       out.addCell(f);
       // Add the vertices
       CCIndex v1 = CCIndexFactory.getIndex();
       out.addCell(v1);
       CCIndex v2 = CCIndexFactory.getIndex();
       out.addCell(v2);

       //CCIndex v3 = CCIndexFactory.getIndex();
       ///out.addCell(v3);
       //CCIndex v4 = CCIndexFactory.getIndex();
       //out.addCell(v4);

       //CCIndex v5 = CCIndexFactory.getIndex();
       //out.addCell(v5);
       //CCIndex v6 = CCIndexFactory.getIndex();
       //out.addCell(v6);

       // I plot only the max and mid anisotropy axis

       CCIndexData V1 = (*indexAttr)[v1];
       CCIndexData V2 = (*indexAttr)[v2];
       //CCIndexData V3 = (*indexAttr)[v3];
       //CCIndexData V4 = (*indexAttr)[v4];
       //CCIndexData V5 = (*indexAttr)[v5];
       //CCIndexData V6 = (*indexAttr)[v6];

       Matrix3d transposeSkewSymm = transpose((*shapeAttr)[f].skewSymmetricTensor);

       V1.pos = V.pos - (0.5 * transposeSkewSymm[0] * AnisotropyVecSize);// * transposeSkewSymm[2][0]);///(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       V2.pos = V.pos + (0.5 * transposeSkewSymm[0] * AnisotropyVecSize);// * transposeSkewSymm[2][0]);///(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
          
       //V3.pos = V.pos - (0.5 * transposeSkewSymm[1] * AnisotropyVecSize * transposeSkewSymm[2][1]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       //V4.pos = V.pos + (0.5 * transposeSkewSymm[1] * AnisotropyVecSize * transposeSkewSymm[2][1]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));

    
       //Point3d minAnisoAxis = transposeSkewSymm[0]^transposeSkewSymm[1];
       //V5.pos = V.pos - (0.5 * minAnisoAxis * AnisotropyVecSize * transposeSkewSymm[2][2]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       //V6.pos = V.pos + (0.5 * minAnisoAxis * AnisotropyVecSize * transposeSkewSymm[2][2]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       

       (*indexAttr)[v1] = V1;
       (*indexAttr)[v2] = V2;
       //(*indexAttr)[v3] = V3;
       //(*indexAttr)[v4] = V4;
       //(*indexAttr)[v5] = V5;
       //(*indexAttr)[v6] = V6;

       //double minAxisNorm = norm(V4.pos - V3.pos);
       //mdxInfo << "minAnisoNorm " << minAxisNorm << endl;
       // Add the edge
       CCIndex e = CCIndexFactory.getIndex();
       out.addCell(e, +v1 -v2);

       //CCIndex e2 = CCIndexFactory.getIndex();
       //out.addCell(e2, +v3 -v4);

       //CCIndex e3 = CCIndexFactory.getIndex();
       //out.addCell(e3, +v5 -v6);

    }
    mesh->drawParms(OutputCC).setGroupVisible("Vertices", true);
    mesh->drawParms(OutputCC).setGroupVisible("Edges", true);
    //mesh->drawParms(OutputCC).setGroupVisible("Faces", true);

    mesh->updateAll(OutputCC);
    return 0;
  }
 
  /*bool ComputeCellShapeQuantifier::initialize(QWidget* parent)
  {
      if(!getProcess(parm("Compute Skew Symmetric tensor process"),anisotropyTensorProcess ))
        throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Compute Skew Symmetric tensor process")));
      if(!getProcess(parm("Compute Antisymmetry tensor process"), antisymmetryTensorProcess))
        throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Compute Antisymmetry tensor process")));
      if(!getProcess(parm("Visualize shape field process"), visualizeCellShapeProcess))
        throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Visualize shape field process")));
      visualizeCellShapeProcess->initialize(parent);
      return true;

  } */
  bool ComputeCellShapeQuantifier::run()
  {
     if(!getProcess(parm("Compute Skew Symmetric tensor process"),anisotropyTensorProcess ))
       throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Compute Skew Symmetric tensor process")));
     if(!getProcess(parm("Compute Antisymmetry tensor process"), antisymmetryTensorProcess))
       throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Compute Antisymmetry tensor process")));
     if(!getProcess(parm("Compute periclinal wall surface ratio for L2 cells"), L2PericlinalSurfRatioProcess))
       throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Compute periclinal wall surface ratio for L2 cells")));
     if(!getProcess(parm("Visualize shape field process"), visualizeCellShapeProcess))
       throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Visualize shape field process")));
     if(!getProcess(parm("Cell volume from heatmap process name"), measureVolumeProcess))
       throw(QString("ComputeCellShapeQuantifier::initialize Cannot make :" + parm("Cell volume from heatmap process name")));

      
     mesh = currentMesh();
     SourceCC = currentMesh()->ccName();
     if(SourceCC.isEmpty())
       throw(QString("No input cell complex specified"));
     if(!mesh->exists(SourceCC))
       throw(QString("Specified input cell complex does not exist"));
     CCStructure &cs = mesh->ccStructure(SourceCC);
     indexAttr = &mesh->indexAttr();
     volumeHeatAttr =  &mesh->heatAttr<double>(parm("Cell Volume Signal Attribute"));
     anisotropyTensorProcess->run();
     antisymmetryTensorProcess->run();
     L2PericlinalSurfRatioProcess->run();
     visualizeCellShapeProcess->run();
     measureVolumeProcess->run(*mesh, cs, *indexAttr, *volumeHeatAttr);
     return true;
  }

  bool WriteCellShapeQuantifier::run()
  {
    mesh = currentMesh();
    QString SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &src = mesh->ccStructure(SourceCC);
    indexAttr = &mesh->indexAttr();
    shapeAttr = &mesh->attributes().attrMap<CCIndex,CellShapeData>(parm("Cell Shape quantifier attribute name"));
    QString volumeHeatName = parm("Cell Volume Signal Attribute");
    //mdxInfo << "Il nome che non trova" << volumeHeatName << endl;
    if(volumeHeatName.isEmpty())
        throw QString("WriteCellShapeQuantifier::run::run Heat map output name is empty");
    volumeHeatAttr =  &mesh->heatAttr<double>(volumeHeatName);
    fileName = parm("File Name");
    
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly))
      throw QString("WriteCellShapeQuantifier::run cannot opened file (%2) for writing").arg(fileName);
    QTextStream out(&file);

    // Write header
    out << QString("Cell Index, Cell Type, Volume, Top Surface Area (top Periclinal cell wall), Bottom Surface Area (bottom Periclinal cell wall), Top/Bottom ratio, Max anisotropy length, Max/Mid anisotropy, Max/Min anisotropy") << endl;
    
    forall(CCIndex f, src.cellsOfDimension(3)) {
      if((*indexAttr)[f].selected)
      {
        out << (*indexAttr)[f].label << ", " << CellTypeToString((*shapeAttr)[f].cellType) << ", " << (*volumeHeatAttr)[(*indexAttr)[f].label] << ", " << (*shapeAttr)[f].topPericlinalWallArea << ", " << (*shapeAttr)[f].bottomPericlinalWallArea << ", " << (*shapeAttr)[f].L2periclinalRatio << ", " << (*shapeAttr)[f].skewSymmetricTensor[0][2]<< ", " << (*shapeAttr)[f].skewSymmetricTensor[0][2]/ (*shapeAttr)[f].skewSymmetricTensor[1][2] << ", " << (*shapeAttr)[f].skewSymmetricTensor[0][2]/ (*shapeAttr)[f].skewSymmetricTensor[2][2] << endl;
      }
    }
    setStatus(QString("Shape analysis CSV file written" ));
  }	  

  REGISTER_PROCESS(FemMembranes);
  REGISTER_PROCESS(FemMembraneSolver);
  REGISTER_PROCESS(FemMembraneRefCfg);
  REGISTER_PROCESS(FemMembraneStressStrain);
  REGISTER_PROCESS(FemMembraneDerivs);
  REGISTER_PROCESS(FemMembraneMaterial);
  REGISTER_PROCESS(FemMembraneCellFacesMaterial);
  REGISTER_PROCESS(FemMembraneAnisoDir);
  REGISTER_PROCESS(FemMembraneAnisoDirMorphogens);
  REGISTER_PROCESS(FemMembranePressure);
  REGISTER_PROCESS(FemMembranePressureDerivs);
  REGISTER_PROCESS(FemMembraneEdgePressure);
  REGISTER_PROCESS(FemMembranePressureEdgeDerivs);
 
  REGISTER_PROCESS(FemMembraneSetGrowth);
  REGISTER_PROCESS(FemMembraneCellFacesGrowth);
  REGISTER_PROCESS(FemMembraneTrichomeProcess);
  REGISTER_PROCESS(FemMembraneSetDirichlet);
  REGISTER_PROCESS(FemMembraneDirichletDerivs);
  REGISTER_PROCESS(FemMembraneGrowth);
  REGISTER_PROCESS(FemMembraneSetGrowthMorphogens);

  REGISTER_PROCESS(FemMembraneBisect);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropyDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropyDirichletDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropySetDirichlet);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropyVisualize);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropySolver);
  REGISTER_PROCESS(FemMembraneCreateAnisotropyMorphogenElement);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthDirichletDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthSetDirichlet);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthVisualize);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthSolver);
  REGISTER_PROCESS(FemMembraneCreateGrowthMorphogenElement);


  REGISTER_PROCESS(FemMembraneVisMaterial);
  REGISTER_PROCESS(FemMembraneVisGrowth);
  REGISTER_PROCESS(FemMembraneComputeCurrentDirections);
  REGISTER_PROCESS(FemMembraneVisDirections);
  //REGISTER_PROCESS(FemMembraneRender);
  REGISTER_PROCESS(FemAnisotropyPropagationFailure);

  REGISTER_PROCESS(AssignCellTypeForShapeQuantifier);

  REGISTER_PROCESS(SkewSymmetricTensor);
  REGISTER_PROCESS(AntiSymmetryTensor); 
  REGISTER_PROCESS(L2PericlinalSurfRatio);
  REGISTER_PROCESS(VisualizeShapeQuantifiers);
  REGISTER_PROCESS(ComputeCellShapeQuantifier);
  REGISTER_PROCESS(WriteCellShapeQuantifier);

}
