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
  class FemMembranes;
  class FemMembraneSolver;
  class FemMembraneGrowth;
  class FemMembraneTrichomeProcess;
  class FemMembraneBisect;

  class L2PericlinalSurfRatio;
  class SkewSymmetricTensor;
  class AntiSymmetryTensor;
  class VisualizeShapeQuantifiers;

  // Main model class
  class FemMembranes : public Process
  {
  public:
    FemMembranes(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/01 FEM Membranes");
      setDesc("FEM simulation with growth and subdivision");

      addParm("Solver Process", "Name of solver process", "Model/CCF/02 FEM Solver");
      addParm("Trichome Process", "Name of the trichome growth assignation process", "Model/CCF/14 Set Trichome Growth");
      addParm("Growth Process", "Name of growth process", "Model/CCF/13 Growth");
      addParm("Subdivide Process", "Name of subdivision process", "Model/CCF/14 Subdivide");
    }

    bool initialize(QWidget *parent);
    bool step();
    bool rewind(QWidget *parent);
    bool finalize(QWidget *parent);
    
    double growthTime = 0;

  private:
    Mesh *mesh = 0;
    QString ccName;
    CCStructure *cs = 0;

    FemMembraneSolver *solverProcess = 0;
    FemMembraneGrowth *growthProcess = 0;
    FemMembraneTrichomeProcess *trichomeProcess = 0;
    FemMembraneBisect *subdivideProcess = 0;

    //double growthTime = 0;
  };

  class FemMembraneSolver : public fem::FemSolver
  {
  public:
    FemMembraneSolver(const Process &proc) : fem::FemSolver(proc) 
    {
      setName("Model/CCF/02 FEM Solver");
      setDesc("FEM Simulation using triangular membrane elements");

      // Update parameters with our own defaults
      setParmDefault("Stress-Strain", "Model/CCF/05 Stress-Strain");

      // Add derivatives processes
      addParm("Element Derivs", "Process for element derivatives", "Model/CCF/03 Triangle Derivs");
      addParm("Pressure Derivs", "Process for pressure derivatives", "Model/CCF/10a Pressure Derivs"); 
      addParm("Pressure Edge Derivs", "Process for pressure derivatives", "Model/CCF/10b Pressure Edge Derivs"); 

      addParm("Dirichlet Derivs", "Process for Dirichlet derivatives", "Model/CCF/16 Dirichlet Derivs");
    }
  };

  class FemMembraneDerivs : public fem::ElementDerivs
  {
  public:
    FemMembraneDerivs(const Process &proc) : ElementDerivs(proc) 
    {
      setName("Model/CCF/03 Triangle Derivs");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneRefCfg : public fem::SetRefCfg
  {
  public:
    FemMembraneRefCfg(const Process &proc) : SetRefCfg(proc) 
    {
      setName("Model/CCF/04 Reference Configuration");

			addParm("Thickness", "Thickness of the membrane element", "1.0");
			setParmDefault("Element Type", "Linear Triangle");
			setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneStressStrain : public fem::StressStrain
  {
  public:
    FemMembraneStressStrain(const Process &proc) : StressStrain(proc) 
    {
      setName("Model/CCF/05 Stress-Strain");

			setParmDefault("Element Type", "Linear Triangle");
			setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneMaterial : public fem::SetTransIsoMaterial
  {
  public:
    FemMembraneMaterial(const Process &proc) : SetTransIsoMaterial(proc) 
    {
      setName("Model/CCF/06 Material Properties");
    }
  };

  
  class FemMembraneCellFacesMaterial : public fem::CellsToFacesTransIsoMaterial
  {
  public:
    FemMembraneCellFacesMaterial(const Process &proc) : CellsToFacesTransIsoMaterial(proc) 
    {
      setName("Model/CCF/07 Cells to Faces Material");
    }
  };

  class FemMembraneAnisoDir : public fem::SetAnisoDir
  {
  public:
    FemMembraneAnisoDir(const Process &proc) : SetAnisoDir(proc) 
    {
      setName("Model/CCF/08a Set Aniso Dir");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneAnisoDirMorphogens : public fem::SetAnisoDirMorphogens
  {
  public:
    FemMembraneAnisoDirMorphogens(const Process &proc) : SetAnisoDirMorphogens(proc) 
    {
      setName("Model/CCF/08b Set Aniso Dir Morphogens");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
      addParm("Morphogen Signal", "Morphogen Signal name",  "E2Morphogen Signal", QStringList() << "E2Morphogen Signal" << "KPArMorphogen Signal" );
    }
  };

  class FemMembranePressure : public fem::SetPressure
  {
  public:
    FemMembranePressure(const Process &proc) : SetPressure(proc) 
    {
      setName("Model/CCF/09a Set Pressure");
    }
  };

  class FemMembranePressureDerivs : public fem::PressureDerivs
  {
  public:
    FemMembranePressureDerivs(const Process &proc) : PressureDerivs(proc) 
    {
      setName("Model/CCF/10a Pressure Derivs");
    }
  };

  class FemMembraneEdgePressure : public fem::SetEdgePressure
  {
  public:
    FemMembraneEdgePressure(const Process &proc) : SetEdgePressure(proc) 
    {
      setName("Model/CCF/09b Set Edge Pressure");
    }
  };

  class FemMembranePressureEdgeDerivs : public fem::PressureEdgeDerivs
  {
  public:
    FemMembranePressureEdgeDerivs(const Process &proc) : PressureEdgeDerivs(proc) 
    {
      setName("Model/CCF/10b Pressure Edge Derivs");
    }
  };



  class FemMembraneSetGrowth : public fem::SetGrowth
  {
  public:
    FemMembraneSetGrowth(const Process &proc) : SetGrowth(proc) 
    {
      setName("Model/CCF/11a Set Growth");
    }
  };

  class FemMembraneSetGrowthMorphogens : public fem::SetGrowthMorphogensFaces
  {
  public:
    FemMembraneSetGrowthMorphogens(const Process &proc) : SetGrowthMorphogensFaces(proc) 
    {
      setName("Model/CCF/11b Set Growth Morphogens (on Faces)");
    }
  };

  class FemMembraneCellFacesGrowth : public fem::CellsToFacesGrowth
  {
  public:
    FemMembraneCellFacesGrowth(const Process &proc) : CellsToFacesGrowth(proc) 
    {
      setName("Model/CCF/12 Cells to Faces Growth");
    }
  };

  class FemMembraneGrowth : public fem::Grow
  {
  public:
    FemMembraneGrowth(const Process &proc) : Grow(proc) 
    {
      setName("Model/CCF/13 Growth");
    }
  };

  class FemMembraneComputeCurrentDirections : public fem::ComputeCurrentDirections
  {
    public:
      FemMembraneComputeCurrentDirections(const Process &proc) : ComputeCurrentDirections(proc)
      {
        setName("Model/CCF/17 Compute Current Directions") ;
      }
  };

  class FemMembraneVisDirections : public fem::VisDirections
  {
    public:
      FemMembraneVisDirections(const Process &proc) : VisDirections(proc)
      {
        setName("Model/CCF/18 Visualize Directions") ;
      }
  };

  //specific trichome growth assignation process
  class FemMembraneTrichomeProcess : public Process
  {
    public:
      FemMembraneTrichomeProcess(const Process &process) : Process(process) 
      {
        setName("Set Trichome Growth");
        setName("Model/CCF/14 Set Trichome Growth");

        setDesc("Set Growth attributes specific for trichome growth. Ovverrides growthIso attributed from previous processes");
  
        addParm("Initial sources Set","Saved set containing tip vertices","InitialSet.txt");
        addParm("Secondary sources Set", "Saved set containing secondary tip vertices", "SecondarySet.txt");
	//addParm("Distance ring-tip", "Distance between the tricome tip and the ring of secondary growth", "20.");
        addParm("Distance boundary-tip", "Distance between the tricome tip and the boundary to be reached before initiationg secondary growth", "40");
        //addParm("Tolerance ring selection", "Tolerance to accept a node into the ring and as a source", "1");
        addParm("Effective distance for growth", "Scaling factor for distance based-growth field", "10.");
        addParm("Effective distance for Secondary Growth", "Effective distance for growth after branching", "5");
        addParm("Decaying factor for distance Secondary Growth", "Decaying factor for the range of action of secondary growth as: EffDistSecondary * e^(-t * Decaying)", "1");
        //addParm("Angular distance", "Angular distance for secondary sources selection (DEG)", "120.");
        addParm("Growth scaling factor", "Global growth rescaling factor", "0.01");
        addParm("FemMembranes Process name", "Name of the main FemMembranes process", "Model/CCF/01 FEM Membranes");
       	addParm("Hodge Attribute",
                "Name of double-valued attribute containing computed Hodge values","Hodge Values");
        addParm("Growth Attribute", "Name of the attribute that holds Growth", "Fem Growth");
        addParm("Distance field 1 Signal", "Name of the signal to visualize distance signal primary sources", "Distance Primary");
        addParm("Distance field 2 Signal", "Name of the signal to visualize distance signal secondary sources", "Distance Secondary");



      }

      /// Set the reference configuration for selected elements
      bool initialize(QWidget *parent);
      bool run();

      FemMembranes *femMembranes = 0;

      private:
       Mesh *mesh = 0;
       QString ccName;
       CCStructure *cs = 0;
       //FemMembranes *femMembranes = 0;

       bool newSources = false;
  };
  class FemMembraneSetDirichlet : public fem::SetDirichlet
  {
  public:
    FemMembraneSetDirichlet(const Process &proc) : SetDirichlet(proc) 
    {
      setName("Model/CCF/15 Set Dirichlet");
    }
  };

  class FemMembraneDirichletDerivs : public fem::DirichletDerivs
  {
  public:
    FemMembraneDirichletDerivs(const Process &proc) : DirichletDerivs(proc) 
    {
      setName("Model/CCF/16 Dirichlet Derivs");
    }
  };
 
  class FemMembraneDiffusionAnisotropySolver : public fem::FemMorphogenSolver
  {
    public:
    FemMembraneDiffusionAnisotropySolver(const Process &proc) :FemMorphogenSolver(proc) 
    {
      setName("Model/CCF/Fem Diffusion Anisotropy/00 Diffusion Solver");

      addParm("Element Diffusion Derivs", "Process for element derivatives", "Model/CCF/Fem Diffusion Anisotropy/01 Diffusion Derivatives");
      addParm("Dirichlet Derivs", "Process for fixed concentration derivatives", "Model/CCF/Fem Diffusion Anisotropy/02 Dirichlet Derivatives");
      setParmDefault("Visualize Process", "Model/CCF/Fem Diffusion Anisotropy/05 Morphogen Visualize");
    }
  };
  class FemMembraneDiffusionAnisotropyDerivs: public fem::MorphogenDiffusionDerivs 
  {
  public:
    FemMembraneDiffusionAnisotropyDerivs(const Process &proc) : MorphogenDiffusionDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Anisotropy/01 Diffusion Derivatives");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");
    }
  };
 
  class FemMembraneDiffusionAnisotropyDirichletDerivs: public fem::MorphogenDirichletDerivs 
  {
  public:
    FemMembraneDiffusionAnisotropyDirichletDerivs(const Process &proc) : MorphogenDirichletDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Anisotropy/02 Dirichlet Derivatives");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "AnisotropyE2Dirichlet", QStringList() << "AnisotropyE2Dirichlet" << "AnisotropyKParDirichlet");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");

    }
  };

  class FemMembraneDiffusionAnisotropySetDirichlet: public fem::SetMorphogenDirichlet 
  {
  public:
    FemMembraneDiffusionAnisotropySetDirichlet(const Process &proc) : SetMorphogenDirichlet(proc) 
    {
       setName("Model/CCF/Fem Diffusion Anisotropy/03 Set Diffusion Dirchlet");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "AnisotropyE2Dirichlet", QStringList() << "AnisotropyE2Dirichlet" << "AnisotropyKParDirichlet");

    }
  };

  class FemMembraneCreateAnisotropyMorphogenElement : public fem::CreateMorphogenElement
  {
     public:
       FemMembraneCreateAnisotropyMorphogenElement(const Process &proc) : CreateMorphogenElement(proc)
       {
          setName("Model/CCF/Fem Diffusion Anisotropy/04 Create Diffusion Element");
       }
  };

  class FemMembraneDiffusionAnisotropyVisualize : public fem::MorphogenVisualize
  {
     public:
       FemMembraneDiffusionAnisotropyVisualize(const Process &proc) : MorphogenVisualize(proc)
       {
          setName("Model/CCF/Fem Diffusion Anisotropy/05 Morphogen Visualize");
          addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");
          addParm("Morphogen Signal", "Morphogen Signal name", "E2Morphogen Signal", QStringList() << "E2Morphogen Signal" << "KParMorphogen Signal" );

       }
  };


////////////////////

  class FemMembraneDiffusionGrowthSolver : public fem::FemMorphogenSolver
  {
    public:
    FemMembraneDiffusionGrowthSolver(const Process &proc) :FemMorphogenSolver(proc) 
    {
      setName("Model/CCF/Fem Diffusion Growth/00 Diffusion Solver");

      addParm("Element Diffusion Derivs", "Process for element derivatives", "Model/CCF/Fem Diffusion Growth/01 Diffusion Derivatives");
      addParm("Dirichlet Derivs", "Process for fixed concentration derivatives", "Model/CCF/Fem Diffusion Growth/02 Dirichlet Derivatives");
      setParmDefault("Visualize Process", "Model/CCF/Fem Diffusion Growth/05 Morphogen Visualize");
    }
  };


  class FemMembraneDiffusionGrowthDerivs: public fem::MorphogenDiffusionDerivs 
  {
  public:
    FemMembraneDiffusionGrowthDerivs(const Process &proc) : MorphogenDiffusionDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Growth/01 Diffusion Derivatives");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");
    }
  };
 
  class FemMembraneDiffusionGrowthDirichletDerivs: public fem::MorphogenDirichletDerivs 
  {
  public:
    FemMembraneDiffusionGrowthDirichletDerivs(const Process &proc) : MorphogenDirichletDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Growth/02 Dirichlet Derivatives");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "GrowthDirichlet1", QStringList() << "GrowthDirichlet1" << "GrowthDirichlet2");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");

    }
  };

  class FemMembraneDiffusionGrowthSetDirichlet: public fem::SetMorphogenDirichlet 
  {
  public:
    FemMembraneDiffusionGrowthSetDirichlet(const Process &proc) : SetMorphogenDirichlet(proc) 
    {
       setName("Model/CCF/Fem Diffusion Growth/03 Set Diffusion Dirchlet");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "GrowthDirichlet1", QStringList() << "GrowthDirichlet1" << "GrowthDirichlet2");

    }
  };

  class FemMembraneCreateGrowthMorphogenElement : public fem::CreateMorphogenElement
  {
     public:
       FemMembraneCreateGrowthMorphogenElement(const Process &proc) : CreateMorphogenElement(proc)
       {
          setName("Model/CCF/Fem Diffusion Growth/04 Create Diffusion Element");
       }
  };

  class FemMembraneDiffusionGrowthVisualize : public fem::MorphogenVisualize
  {
     public:
       FemMembraneDiffusionGrowthVisualize(const Process &proc) : MorphogenVisualize(proc)
       {
          setName("Model/CCF/Fem Diffusion Growth/05 Morphogen Visualize");
          addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");
          addParm("Morphogen Signal", "Morphogen Signal name", "Growth1 Morph Signal", QStringList() << "Growth1 Morph Signal" << "Growth2 Morph Signal" );

       }
  };



  /// Subdivide object
  class FemMembraneSubdivide : public mdx::Subdivide
  {
  public:
    FemMembraneSubdivide() {}

    FemMembraneSubdivide(Mesh &mesh, fem::ElasticTriangle3Attr &elementAttr, fem::Triangle3Attr<Point1d> &elementMorphoAnisoAttr, fem::Triangle3Attr<Point1d> &elementMorphoGrowthAttr,  
        fem::TransIsoMaterialAttr &materialAttr, fem::PressureAttr &pressureAttr, fem::PressureEdgeAttr &pressureEdgeAttr, fem::DirichletAttr &dirichletAttr, fem::GrowthAttr &growthAttr, fem::MorphogenDataAttr &morphogenE2DataAttr,
        fem::MorphogenDataAttr &morphogenKParDataAttr,  
        fem::MorphogenDataAttr &morphogenGrowth1DataAttr, fem::MorphogenDataAttr &morphogenGrowth2DataAttr, fem::MorphogenDirichletAttr &moprhogenDirichletE2Attr, 
        fem::MorphogenDirichletAttr &moprhogenDirichletKParAttr, fem::MorphogenDirichletAttr &moprhogenDirichletGrowth1Attr, fem::MorphogenDirichletAttr &moprhogenDirichletGrowth2Attr)  : 
           mdxSubdivide(mesh), elementSubdivide(mesh.indexAttr(), elementAttr), elementMorphogAnisoSubdivide(mesh.indexAttr(), elementMorphoAnisoAttr),  elementMorphogGrowthSubdivide(mesh.indexAttr(), elementMorphoGrowthAttr), 
        materialSubdivide(materialAttr), pressureSubdivide(pressureAttr), pressureEdgeSubdivide(pressureEdgeAttr), dirichletSubdivide(dirichletAttr), growthSubdivide(growthAttr), morphogenE2Subdivide(mesh.indexAttr(), morphogenE2DataAttr), morphogenKParSubdivide(mesh.indexAttr(), morphogenKParDataAttr),
        morphogenGrowth1Subdivide(mesh.indexAttr(), morphogenGrowth1DataAttr), morphogenGrowth2Subdivide(mesh.indexAttr(), morphogenGrowth2DataAttr), morphDirE2Subdivide(mesh.indexAttr(), moprhogenDirichletE2Attr), morphDirKParSubdivide(mesh.indexAttr(), moprhogenDirichletKParAttr), morphDirGrowth1Subdivide(mesh.indexAttr(), moprhogenDirichletGrowth1Attr), morphDirGrowth2Subdivide(mesh.indexAttr(), moprhogenDirichletGrowth2Attr) {}

    // Method to split the element data
    void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct &ss, 
        CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(), double interpPos = 0.5) 
    {
      mdxSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

      // Propagate the material parameters
      if(dim == 2 or dim == 3) {
        elementSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        materialSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        pressureSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        pressureEdgeSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

        growthSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

        elementMorphogAnisoSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        elementMorphogGrowthSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

      
        morphogenE2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphogenKParSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphogenGrowth1Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphogenGrowth2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        
      }
      else if(dim == 1)
      {
        dirichletSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphDirKParSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphDirE2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

        morphDirGrowth1Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphDirGrowth2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
       
      }
    }
    MDXSubdivide mdxSubdivide;
    fem::ElasticTriangle3::Subdivide elementSubdivide;
    fem::Triangle3<Point1d>::Subdivide elementMorphogAnisoSubdivide;
    fem::Triangle3<Point1d>::Subdivide elementMorphogGrowthSubdivide;

    fem::TransIsoMaterial::Subdivide materialSubdivide;
    fem::Pressure::Subdivide pressureSubdivide;
    fem::PressureEdge::Subdivide pressureEdgeSubdivide;
    fem::Growth::Subdivide growthSubdivide;
    fem::Dirichlet::Subdivide dirichletSubdivide;

    fem::MorphogenData::Subdivide morphogenE2Subdivide;
    fem::MorphogenData::Subdivide morphogenKParSubdivide;
    fem::MorphogenDirichlet::Subdivide morphDirE2Subdivide;
    fem::MorphogenDirichlet::Subdivide morphDirKParSubdivide;

    fem::MorphogenData::Subdivide morphogenGrowth1Subdivide;
    fem::MorphogenDirichlet::Subdivide morphDirGrowth1Subdivide;
    fem::MorphogenData::Subdivide morphogenGrowth2Subdivide;
    fem::MorphogenDirichlet::Subdivide morphDirGrowth2Subdivide;


  };

  class FemMembraneBisect : public SubdivideBisectTriangle
  {
  public:
    FemMembraneBisect(const Process &proc) : SubdivideBisectTriangle(proc) 
    {
      setName("Model/CCF/14 Subdivide");

      addParm("Element Attribute", "Attribute to store nodal values", "Triangle Element");
      
      addParm("Material Attribute", "Name of the attribute that holds material properties", "TransIso Material");
      addParm("Pressure Attribute", "Name of the attribute that holds pressure", "Fem Pressure");
      addParm("Pressure Edge Attribute", "Name of the attribute that holds pressure", "Fem Edge Pressure");

      addParm("Growth Attribute", "Name of the attribute that holds growth", "Fem Growth");
      addParm("Dirichlet Attribute", "Name of the attribute that holds dirichlet", "Fem Dirichlet");
      addParm("Morphogen Anisotropy Element Attribute", "Attribute to store morphogenetic element for anisotropy", "Morphogen Triangle Element Growth");
      addParm("Morphogen Growth Element Attribute", "Attribute to store morphogenetic element for growth", "Morphogen Triangle Element Growth");

      addParm("Morphogen Anisotropy E2 Data", "Attribute to store morphogenetic concentration for E2", "AnisotropyE2Morphogen");
      addParm("Morphogen Anisotropy KPar Data", "Attribute to store morphogenetic concentration for KPar", "AnisotropyKParMorphogen");
      addParm("Morphogen Anisotropy E2 Dirichlet", "Attribute to store morphogenetic Dirichlet for E2", "AnisotropyE2Dirichlet");
      addParm("Morphogen Anisotropy KPar Dirichlet", "Attribute to store morphogenetic Dirichlet for KPar", "AnisotropyKParDirichelt");

      addParm("Morphogen Growth Data 1", "Attribute to store morphogenetic concentration for growth1", "GrowthMorphogen1");
      addParm("Morphogen Growth Data 2", "Attribute to store morphogenetic concentration for growth2", "GrowthMorphogen2");
      addParm("Morphogen Growth Dirichlet 1", "Attribute to store morphogenetic Dirichlet for growth1", "GrowthDirichlet1");
      addParm("Morphogen Growth Dirichlet 2", "Attribute to store morphogenetic Dirichlet for growth2", "GrowthDirichlet2");
  
    }
    using SubdivideBisectTriangle::run;
    bool run();
  };

  class FemMembraneVisMaterial : public fem::VisTransIsoMaterial
  {
  public:
    FemMembraneVisMaterial(const Process &proc) : VisTransIsoMaterial(proc) 
    {
      setName("Model/CCF/20 Visualize Material");
    }
  };

  class FemMembraneVisGrowth : public fem::VisGrowth
  {
  public:
    FemMembraneVisGrowth(const Process &proc) : VisGrowth(proc) 
    {
      setName("Model/CCF/21 Visualize Growth");
    }
  };


  //class FemMembraneRender : public fem::VectorRender
  //{
  //public:
  //  FemMembraneRender(const Process &proc) : VectorRender(proc) 
  //  {
  //    setName("Model/CCF/22 Vector Render");
  //  }
  //};
  class FemAnisotropyPropagationFailure : public fem::DisplayFailedAnisotropyPropagation
  { 
    public:
    FemAnisotropyPropagationFailure(const Process &proc) : DisplayFailedAnisotropyPropagation(proc) 
    {
      setName("Model/CCF/22 Display Anisotropy Propagation Failure");
    }
  };
 
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

        addParm("Compute Skew Symmetric tensor process", "Compute Skew Symmetric tensor process", "Model/CCF/60 Shape Quantifier/01 Compute Skew Symmetric Tensor");
        addParm("Compute Antisymmetry tensor process", "Compute Antisymmetry tensor process", "Model/CCF/60 Shape Quantifier/02 Compute Antisymmetry Tensor");
        addParm("Compute periclinal wall surface ratio for L2 cells", "For each L2 cell, computes the ratio of surface area shared with top L1 cells and the surface area shared with L3 layer",  "Model/CCF/60 Shape Quantifier/03 Compute L2 periclinal surface ratio for selected ovule zone");
        addParm("Cell volume from heatmap process name", "Compute cell volume from heatmap process, provide the process name", "Mesh/Heat Map/Measures3D/Geometry/Volume");
        addParm("Cell Volume Signal Attribute", "Cell Volume Signal Attribute as from Volume heatmap computation, Mesh/Heat Map/Measures3D/Geometry/Volume", "Volume");
        addParm("Visualize shape field process", "Visualize shape field", "Model/CCF/60 Shape Quantifier/04 Visualize Shape Field");
        //addParm("Write Periclinal wall ratio to a file for selected cells", "Writes periclinal wall ratio, top wall size, bottom wall size, general cell label, cell label, cell volume", "testShapeQuantifier.csv");

      }
      //bool initialize(QWidget* parent);
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      IntDoubleAttr *volumeHeatAttr = 0;
      //AuxinGradient::CellDataAttr *cellAttr = 0;
      //AuxinGradient::EdgeDataAttr *edgeAttr = 0;
      //SkewSymmetricTensor *skewSymmetricTensorProcess = 0;
      //CellShapeAttr *shapeAttr = 0;
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
        setName("Model/CCF/60 Shape Quantifier/01 Compute Skew Symmetric Tensor");

        setDesc("Compute cell shape anisotropy.");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      

        //addParm("Tissue Process", "Name of process for Cell Tissue simulation", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");


      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      //CellTissueProcess *tissueProcess = 0;

      //AuxinGradient::CellDataAttr *cellAttr = 0;
      //AuxinGradient::EdgeDataAttr *edgeAttr = 0;

    private:

  };

  class AntiSymmetryTensor : public Process
  {
    public:
      AntiSymmetryTensor(const Process &process) : Process(process) 
      {
        setName("Model/CCF/60 Shape Quantifier/02 Compute Antisymmetry Tensor");

        setDesc("Compute cell shape anisotropy.");
        setIcon(QIcon(":/images/CellType.png"));

        addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      //AuxinGradient::CellDataAttr *cellAttr = 0;
      //AuxinGradient::EdgeDataAttr *edgeAttr = 0;

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
        //addParm("Select L1 dome cells", "Select L1 dome cells and label as L1Dome", "L1 Dome");

      }
      bool run();
      CCIndexDataAttr  *indexAttr = 0;
      CellShapeAttr *shapeAttr = 0;
      Mesh *mesh = 0;
      QString SourceCC;
      double spuriousContactCutoff;
      //AuxinGradient::CellDataAttr *cellAttr = 0;
      //AuxinGradient::EdgeDataAttr *edgeAttr = 0;

    private:

  };

  class VisualizeShapeQuantifiers : public Process
  {
    public:
    //enum ParmNames { pOutputCC, pVectorSize, pNumParms };

    VisualizeShapeQuantifiers(const Process &process) : Process(process) 
    {
      setName("Model/CCF/60 Shape Quantifier/04 Visualize Shape Field");
      setDesc("Draw shape fields and anisotropy vectors.");
      setIcon(QIcon(":/images/Default.png"));


      addParm("Shape Attribute Name", "Shape Attribute Name", "CellShapeData");
      //addParm("Source CC", "Name of source cell complex", "");
      //addParm("Visualized field", "Which shape field to visualize", "Anisotropy", QStringList() << "Anisotropy" << "Max Asymmetry" << "Min Asymmetry" << "Anisotropy + Max Asymmetry" << "L2 Periclinal Wall Ratio" <<  "None" );
      addParm("Output CC", "Name of output cell complex", "Draw cell anisotropy");
      addParm("Anisotropy Vector Size", "Amount to scale anisotropy vector", "1.0");
    }
    //bool initialize(QWidget* parent);
  
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
    //Mesh *mesh = 0;

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
      //bool initialize(QWidget* parent);
      bool run();
      Mesh *mesh = 0;
      CCIndexDataAttr  *indexAttr = 0;
      //AuxinGradient::CellDataAttr *cellAttr = 0;
      //AuxinGradient::EdgeDataAttr *edgeAttr = 0;
      //SkewSymmetricTensor *skewSymmetricTensorProcess = 0;
      CellShapeAttr *shapeAttr = 0;
      IntDoubleAttr *volumeHeatAttr = 0;
      QString fileName; 
    private:

  };


  
}
#endif

