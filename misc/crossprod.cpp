using Eigen::Upper;
typedef Eigen::Map<Eigen::MatrixXd>  MMat;

MMat V(as<MMat> (RV));
MMat VtV(as<MMat> (RVtV));

VtV.setZero().selfadjointView<Upper>().rankUpdate(V.adjoint());

