#include<algorithm>
#include<iostream>
#include<vector>
#include<tuple>
#include<utility>
#include<array>



template<typename Target,typename Tuple, int N, bool end>
struct Extractor {
  template < typename ... Args >
    static Target index(int i,const Tuple& t, Args && ... args)
    {
      return Extractor<Target, Tuple,  N+1, 
             std::tuple_size<Tuple>::value == N+1>::
               index(i, t , std::forward<Args>(args)..., std::get<N>(t).vec[i] );
    }
};

template < typename Target, typename Tuple, int N >
struct Extractor<Target,Tuple,N,true>
{
  template < typename ... Args >
    static Target index(int i,Tuple const& t, 
        Args && ... args) { 
      return Target(std::forward<Args>(args)...); }
};

template < typename ... Vs > 
  std::tuple<typename std::remove_reference<Vs>::type::type_t...>
extract(int i , const std::tuple<Vs...>& tpl)
{
  return Extractor<std::tuple<typename std::remove_reference<Vs>::type::type_t...>,
         std::tuple<Vs...>, 0,
         std::tuple_size<std::tuple<Vs...> >::value == 0>::index(i,tpl);
}


struct Vec {
  std::array<float,5> vec;
  typedef float type_t;

  template<typename... OPs,typename... VALs>
    Vec& operator=(const std::pair< std::tuple<VALs...> , std::tuple<OPs...> >& e) {
      for( int i = 0 ; i < vec.size() ; ++i ) {
        vec[i] = eval( extract(i,e.first) , e.second );
      }
    }
};




template<int OpPos,int ValPos, bool end>
struct StackMachine {
  template<typename... OPs,typename... VALs>
    static void eval_pos( std::tuple<VALs...>& vals , const std::tuple<OPs...> & ops )
    {
      std::get<ValPos+1>( vals ) =
        std::get<OpPos>(ops).apply( std::get<ValPos>( vals ) , 
            std::get<ValPos+1>( vals ) );
      StackMachine<OpPos+1,ValPos+1,sizeof...(OPs) == OpPos+1>::eval_pos(vals,ops);
    }
};

template<int OpPos,int ValPos>
struct StackMachine<OpPos,ValPos,true> {
  template<typename... OPs,typename... VALs>
    static void eval_pos( std::tuple<VALs...>& vals , 
        const std::tuple<OPs...> & ops )
    {}
};



  template<typename... OPs,typename... VALs>
float eval( const std::tuple<VALs...>& vals , const std::tuple<OPs...> & ops )
{
  StackMachine<0,0,false>::eval_pos(const_cast<std::tuple<VALs...>&>(vals),ops);
  return std::get<sizeof...(OPs)>(vals);
}




struct OpMul {
  static int apply(const float& lhs,const float& rhs)  {
    return lhs*rhs;
  }
};

  std::pair< std::tuple< const Vec&, const Vec& > , std::tuple<OpMul> >
operator*(const Vec& lhs,const Vec& rhs)
{
  return std::make_pair( std::tuple< const Vec&, const Vec& >( lhs , rhs ) , 
      std::tuple<OpMul>( OpMul() ) );
}

template<typename... OPs,typename... VALs>
  std::pair< std::tuple< const Vec&, VALs... > , std::tuple<OPs...,OpMul> >
operator*(const Vec& lhs,const std::pair< std::tuple< VALs... > , std::tuple<OPs...> >& rhs)
{
  return std::make_pair( std::tuple_cat( rhs.first , std::tuple< const Vec& >(lhs)  ) , 
      std::tuple_cat( rhs.second , std::tuple<OpMul>( OpMul() )  ) );
}

template<typename... OPs,typename... VALs>
  std::pair< std::tuple< const Vec&, VALs... > , std::tuple<OPs...,OpMul> >
operator*(const std::pair< std::tuple< VALs... > , std::tuple<OPs...> >& lhs,
    const Vec& rhs)
{
  return std::make_pair( std::tuple_cat( lhs.first , std::tuple< const Vec& >(rhs)  ) , 
      std::tuple_cat( lhs.second , std::tuple<OpMul>( OpMul() ) ) );
}

int main(int argc, char * argv[])
{
  Vec d,c,b,a;


  for( int i = 0 ; i < d.vec.size() ; ++i ) {
    a.vec[i] = 10+i+atoi(argv[2]);
    b.vec[i] = 20+i+atoi(argv[5]);
    c.vec[i] = 30+i+atoi(argv[4]);
    d.vec[i] = 0+atoi(argv[3]);
  }

  asm("#test1");
  d = a * b * c * a;
  asm("#test2");

  for( int i = 0 ; i < d.vec.size() ; ++i ) 
    std::cout << d.vec[i] << std::endl;
}

