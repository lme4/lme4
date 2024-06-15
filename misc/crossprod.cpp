using Eigen::MatrixXd;
using Eigen::Upper;
typedef Eigen::Map<MatrixXd>  MMat;

MMat V(as<MMat> (RV));
MMat VtV(as<MMat> (RVtV));

VtV.setZero().selfadjointView<Upper>().rankUpdate(V.adjoint());

