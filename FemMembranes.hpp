#ifndef FEM_MEMBRANES_HPP
#define FEM_MEMBRANES_HPP

#include <MDXProcessFem.hpp>
#include <MDXProcessFemMorphogen.hpp>
#include <MeshProcessSystem.hpp>
#include <MeshProcessStructure.hpp>
#include <MDXProcessTissue.hpp>
#include <MDXProcessCellDivide.hpp>
#include <Attributes.hpp>
#include <MeshProcessSystemRender.hpp>
#include <MeshUtils.hpp>
#include <MeshProcessSelection.hpp>
#include <MeshProcessMeasures3D.hpp>
#include <Triangle3.hpp>

#include <Solver.hpp>
namespace mdx
{

  class L2PericlinalSurfRatio;
  class SkewSymmetricTensor;
  class AntiSymmetryTensor;
  class VisualizeShapeQuantifiers;

   
  enum CellType {L1Dome, L1, L2, L3, pSMC, CC, genericCell};
  static CellType stringToCellType(const QString &str)
    {
      if(str == "L1 Dome")
        return(L1Dome);
      else if(str == "L1")
        return(L1);
      else if(str == "L2")
        return(L2);
      else if(str=="L3")
        return(L3);
      else if(str=="pSMC")
        return(pSMC);
      else if(str=="CC")
        return(CC);
      else if(str== "Generic Cell")
        return(genericCell);
      else 
        throw(QString("Bad cell type %1").arg(str));
    }

  static QString CellTypeToString(const CellType &cellType)
    {
      if(cellType == CellType::L1Dome)
        return("L1 Dome");
      else if(cellType == CellType::L1)
        return("L1");
      else if(cellType == CellType::L2)
        return("L2");
      else if(cellType == CellType::L3)
        return("L3");
      else if(cellType == CellType::pSMC)
        return("pSMC");
      else if(cellType == CellType::CC)
        return("CC");
      else if(cellType == CellType::genericCell)
        return("GenericCell");
      else 
        throw(QString("Bad cell type %1").arg(cellType));
    }

  struct CellShapeData
  {
     Matrix3d skewSymmetricTensor; 
     Point3d asymmetry;
     CellType cellType = genericCell;
     double L2periclinalRatio = -1;
     double bottomPericlinalWallArea = -1;
     double topPericlinalWallArea = -1;
     
     CellShapeData() {}
  
     bool operator==(const CellShapeData &other) const
     {
      if(skewSymmetricTensor == other.skewSymmetricTensor and asymmetry== other.asymmetry and cellType == other.cellType and L2periclinalRatio == other.L2periclinalRatio and bottomPericlinalWallArea == other.bottomPericlinalWallArea and topPericlinalWallArea == other.topPericlinalWallArea)
        return true;
      return false;
     }      
  };
  typedef AttrMap<CCIndex,CellShapeData> CellShapeAttr;

  bool inline readAttr(CellShapeData &m, const QByteArray &ba, size_t &pos) 
  { 
    return mdx::readChar((char *)&m, sizeof(CellShapeData), ba, pos);
  }
  bool inline writeAttr(const CellShapeData  &m, QByteArray &ba)
  { 
    return mdx::writeChar((char *)&m, sizeof(CellShapeData), ba);
  }

  class AssignCellTypeForShapeQuantifier : public Process
  {
    public:
      AssignCellTypeForShapeQuantifier (const Process &process) : Process(process) 
      {
        setName("Model/CCF/50 Set Cell Type");
        setDesc("Assign cell type for shape quantifiers (L1 dome and L1 cells).");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Assign cell type for selected cells", "Assign cell type for selected cells, L1 and L1 Dome are required for shape quantification", "L1", QStringList()<< "L1" << "L2" << "pSMC" << "CC");
        addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");

      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      CellType cellType;
  };

  class ComputeCellShapeQuantifier : public Process
  {
    public:
      ComputeCellShapeQuantifier(const Process &process) : Process(process) 
      {
        setName("Model/CCF/60 Shape Quantifier/00 Global Shape Quantifier Process");
        setDesc("Compute cell shape anisotropy and antisymmetry.");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Compute Skew Symmetric tensor process", "Compute Skew Symmetric tensor process", "Model/CCF/60 Shape Quantifier/01 Compute Covariance Matrix");
        addParm("Compute Antisymmetry tensor process", "Compute Antisymmetry tensor process", "Model/CCF/60 Shape Quantifier/02 Compute Antisymmetry Matrix");
        addParm("Compute periclinal wall surface ratio for L2 cells", "For each L2 cell, computes the ratio of surface area shared with top L1 cells and the surface area shared with L3 layer",  "Model/CCF/60 Shape Quantifier/03 Compute L2 periclinal surface ratio for selected ovule zone");
        addParm("Cell volume from heatmap process name", "Compute cell volume from heatmap process, provide the process name", "Mesh/Heat Map/Measures3D/Geometry/Volume");
        addParm("Cell Volume Signal Attribute", "Cell Volume Signal Attribute as from Volume heatmap computation, Mesh/Heat Map/Measures3D/Geometry/Volume", "Volume");
        addParm("Visualize shape field process", "Visualize shape field", "Model/CCF/60 Shape Quantifier/04 Visualize Shape Field");
        addParm("Write Periclinal wall ratio to a file for selected cells", "Writes periclinal wall ratio, top wall size, bottom wall size, general cell label, cell label, cell volume", "testShapeQuantifier.csv");

      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      IntDoubleAttr *volumeHeatAttr = 0;
      SkewSymmetricTensor *anisotropyTensorProcess = 0;
      AntiSymmetryTensor *antisymmetryTensorProcess = 0;
      L2PericlinalSurfRatio *L2PericlinalSurfRatioProcess = 0;
      VisualizeShapeQuantifiers *visualizeCellShapeProcess = 0;
      MeasureVolume *measureVolumeProcess = 0; 
    private:

  };

  class SkewSymmetricTensor : public Process
  {
    public:
      SkewSymmetricTensor(const Process &process) : Process(process) 
      {
        setName("Model/CCF/60 Shape Quantifier/01 Compute Covariance Matrix");

        setDesc("Compute cell shape anisotropy through covariance matrix.");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      



      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;

    private:

  };

  class AntiSymmetryTensor : public Process
  {
    public:
      AntiSymmetryTensor(const Process &process) : Process(process) 
      {
        setName("Model/CCF/60 Shape Quantifier/02 Compute Antisymmetry Matrix");

        setDesc("Compute cell shape antisymmetry.");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;

    private:

  };


  class L2PericlinalSurfRatio : public Process
  {
    public:
      L2PericlinalSurfRatio(const Process &process) : Process(process) 
      {
        setName("Model/CCF/60 Shape Quantifier/03 Compute L2 periclinal surface ratio for selected ovule zone");
        setDesc("For each L2 cell, computes the ratio of surface area shared with L1 cells and the surface area shared with L3 layer");
       
        addParm("Cutoff value for spurious contact", "Set lower cutoff value to filter out spurious contacts as just edge contact or point contact", "1.");
        addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      

        setIcon(QIcon(":/images/CellType.png"));

      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      double spuriousContactCutoff;

    private:

  };

  class VisualizeShapeQuantifiers : public Process
  {
    public:

    VisualizeShapeQuantifiers(const Process &process) : Process(process) 
    {
      setName("Model/CCF/60 Shape Quantifier/04 Visualize Shape Field");
      setDesc("Draw shape fields and anisotropy vectors.");
      setIcon(QIcon(":/images/Default.png"));


      addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      addParm("Output CC", "Name of output cell complex", "Draw cell anisotropy");
      addParm("Anisotropy Vector Size", "Amount to scale anisotropy vector", "1.0");
      addParm("Visualize Max Anisotropy Vector", "Visualize Max Anisotropy Vector", "True", QStringList() << "True" << "False");
      addParm("Visualize Mid Anisotropy Vector", "Visualize Mid Anisotropy Vector", "True", QStringList() << "True" << "False");
      addParm("Visualize Min Anisotropy Vector", "Visualize Min Anisotropy Vector", "True", QStringList() << "True" << "False");

    }
  
    bool run() { return run(currentMesh()); }
    bool run(Mesh *mesh);

  private:
    // Parameters
    QString SourceCC;
    QString OutputCC;
    bool DrawAnisoVec;
    double AnisotropyVecSize;
    CellShapeAttr *shapeAttr = 0;
    CCIndexDataAttr *indexAttr = 0;

  };

  class WriteCellShapeQuantifier : public Process
  {
    public:
      WriteCellShapeQuantifier(const Process &process) : Process(process) 
      {
        setName("Model/CCF/70 Write Shape Quantifier");
        setDesc("Write periclinal wall ratio, bottom wall surface, top wall surface, cell volume, for selected L2 cells with their labeling and cell typ (L2 or pSMC) ");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Cell Shape quantifier attribute name", "Cell Shape quantifier attribute name", "CellShapeData");
        addParm("Cell Volume Signal Attribute", "Cell Volume Signal Attribute as from Volume heatmap computation, Mesh/Heat Map/Measures3D/Geometry/Volume", "Volume");
        addParm("File Name", "File name to write data", "Ovule_EM_C_Test.csv");
       

      }
      bool run();
      Mesh *mesh = 0;
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      IntDoubleAttr *volumeHeatAttr = 0;
      QString fileName; 
    private:

  };


  
}
#endif

