#include "FemMembranes.hpp"

namespace mdx
{

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
    CCIndexVec volumes = selectedVolumes(cs, *indexAttr);
    if (volumes.size() == 0) {
      forall(CCIndex f, cs.cellsOfDimension(3))
        volumes.push_back(f);
    }
    for(uint i = 0; i < volumes.size(); i++) {
       CCIndex f = volumes[i];
       if ((*indexAttr)[f].selected == true){   
        (*shapeAttr)[f].cellType = cellType;
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
      (*shapeAttr)[c].skewSymmetricTensor =  transpose(eigVect);
      (*shapeAttr)[c].skewSymmetricTensor[0][2] = eigVal[0];
      (*shapeAttr)[c].skewSymmetricTensor[1][2] = eigVal[1];
      (*shapeAttr)[c].skewSymmetricTensor[2][2] = eigVal[2];
        
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
    //create a vector of L1 and L1 Dome cells
    // and a vector of L2 cells, to save computational time
    std::vector <CCIndex> L1Cells;
    //std::vector <CCIndex> L1DomeCells;
    std::vector <CCIndex> L2Cells;

    for(CCIndex c : cs.cellsOfDimension(3)){
       if((*shapeAttr)[c].cellType == 1)
         L1Cells.push_back(c);
  
    }

    if (L1Cells.size() == 0 /*or L1DomeCells.size() == 0*/)
      throw(QString("60 Shape Quantifier/03 Compute L2 periclinal surface ratio::run no L1  cells selected. Plese run the process: 50 Set Cell Type and assign L1 cells"));
    
    //assign from L1 and L1 Dome neighborhood relation, L2 cells (only if they belong to the selected area).
    //also assign the L2 attribute to those cells
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
    CCIndexVec volumes = selectedVolumes(src, *indexAttr);
    if (volumes.size() == 0) {
      forall(CCIndex f, src.cellsOfDimension(3))
        volumes.push_back(f);
    }
    for(uint i = 0; i < volumes.size(); i++) {
       CCIndex f = volumes[i];

       //as a scalar I plot the max Dimension/ mid Dimension
       (anisoAttr)[f] = ((*shapeAttr)[f].skewSymmetricTensor[0][2] / ((*shapeAttr)[f].skewSymmetricTensor[0][2] + (*shapeAttr)[f].skewSymmetricTensor[1][2] + (*shapeAttr)[f].skewSymmetricTensor[2][2]));  
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
       //out.addCell(v1);
       CCIndex v2 = CCIndexFactory.getIndex();
       //out.addCell(v2);

       CCIndex v3 = CCIndexFactory.getIndex();
       //out.addCell(v3);
       CCIndex v4 = CCIndexFactory.getIndex();
       //out.addCell(v4);

       CCIndex v5 = CCIndexFactory.getIndex();
       //out.addCell(v5);
       CCIndex v6 = CCIndexFactory.getIndex();
       //out.addCell(v6);

       // I plot only the max and mid anisotropy axis

       CCIndexData V1 = (*indexAttr)[v1];
       CCIndexData V2 = (*indexAttr)[v2];
       CCIndexData V3 = (*indexAttr)[v3];
       CCIndexData V4 = (*indexAttr)[v4];
       CCIndexData V5 = (*indexAttr)[v5];
       CCIndexData V6 = (*indexAttr)[v6];

       Matrix3d transposeSkewSymm = transpose((*shapeAttr)[f].skewSymmetricTensor);

       V1.pos = V.pos - (0.5 * transposeSkewSymm[0] * AnisotropyVecSize * transposeSkewSymm[2][0]/(transposeSkewSymm[2][0]  + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       V2.pos = V.pos + (0.5 * transposeSkewSymm[0] * AnisotropyVecSize * transposeSkewSymm[2][0]/(transposeSkewSymm[2][0]  + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
          
       V3.pos = V.pos - (0.5 * transposeSkewSymm[1] * AnisotropyVecSize * transposeSkewSymm[2][1]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       V4.pos = V.pos + (0.5 * transposeSkewSymm[1] * AnisotropyVecSize * transposeSkewSymm[2][1]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));

    
       Point3d minAnisoAxis = transposeSkewSymm[0]^transposeSkewSymm[1];
       V5.pos = V.pos - (0.5 * minAnisoAxis * AnisotropyVecSize * transposeSkewSymm[2][2]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       V6.pos = V.pos + (0.5 * minAnisoAxis * AnisotropyVecSize * transposeSkewSymm[2][2]/(transposeSkewSymm[2][0] + transposeSkewSymm[2][1] + transposeSkewSymm[2][2]));
       

       (*indexAttr)[v1] = V1;
       (*indexAttr)[v2] = V2;
       (*indexAttr)[v3] = V3;
       (*indexAttr)[v4] = V4;
       (*indexAttr)[v5] = V5;
       (*indexAttr)[v6] = V6;

       //double minAxisNorm = norm(V4.pos - V3.pos);
       //mdxInfo << "minAnisoNorm " << minAxisNorm << endl;
       // Add the edge

       CCIndex e = CCIndexFactory.getIndex();
       CCIndex e2 = CCIndexFactory.getIndex();
       CCIndex e3 = CCIndexFactory.getIndex();

       if(parm("Visualize Max Anisotropy Vector") == "True"){
         out.addCell(v1);
         out.addCell(v2);
         out.addCell(e, +v1 -v2);
       }
       if(parm("Visualize Mid Anisotropy Vector") == "True"){
         out.addCell(v3);
         out.addCell(v4);
         out.addCell(e2, +v3 -v4);
       }
       if(parm("Visualize Min Anisotropy Vector") == "True"){
         out.addCell(v5);
         out.addCell(v6);
         out.addCell(e3, +v5 -v6);
       }
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
    out << QString("Cell Index, Cell Type, Volume, Top Surface Area (top Periclinal cell wall), Bottom Surface Area (bottom Periclinal cell wall), Top/Bottom ratio, Max anisotropy length, Max/Mid anisotropy, Max/Min anisotropy, Max anisotropy eigenVec x, Max eigenVec anisotropy y, Max eigenVec anisotropy z") << endl;
    
    forall(CCIndex f, src.cellsOfDimension(3)) {
      if((*indexAttr)[f].selected)
      {
        out << (*indexAttr)[f].label << ", " << CellTypeToString((*shapeAttr)[f].cellType) << ", " << (*volumeHeatAttr)[(*indexAttr)[f].label] << ", " << (*shapeAttr)[f].topPericlinalWallArea << ", " << (*shapeAttr)[f].bottomPericlinalWallArea << ", " << (*shapeAttr)[f].L2periclinalRatio << ", " << (*shapeAttr)[f].skewSymmetricTensor[0][2]<< ", " << (*shapeAttr)[f].skewSymmetricTensor[0][2]/ (*shapeAttr)[f].skewSymmetricTensor[1][2] << ", " << (*shapeAttr)[f].skewSymmetricTensor[0][2]/ (*shapeAttr)[f].skewSymmetricTensor[2][2] << ", " << (*shapeAttr)[f].skewSymmetricTensor[0][0] <<  ", " << (*shapeAttr)[f].skewSymmetricTensor[1][0] <<  ", " << (*shapeAttr)[f].skewSymmetricTensor[2][0] <<endl;
      }
    }
    setStatus(QString("Shape analysis CSV file written" ));
  }	  


  REGISTER_PROCESS(AssignCellTypeForShapeQuantifier);

  REGISTER_PROCESS(SkewSymmetricTensor);
  REGISTER_PROCESS(AntiSymmetryTensor); 
  REGISTER_PROCESS(L2PericlinalSurfRatio);
  REGISTER_PROCESS(VisualizeShapeQuantifiers);
  REGISTER_PROCESS(ComputeCellShapeQuantifier);
  REGISTER_PROCESS(WriteCellShapeQuantifier);

}
