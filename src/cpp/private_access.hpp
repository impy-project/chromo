//************************************************************************************
//
// This file is a modified version of cpp-member-accessor
// Project home: https://github.com/hliberacki/cpp-member-accessor
// This header is released under the MIT license, see license text below.
//
// Copyright Hubert Liberacki (hliberacki@gmail.com)
// Copyright Krzysztof Ostrowski
//
// MIT License
//
// Copyright (c) 2018 Hubert Liberacki
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//************************************************************************************

#ifndef PRIVATE_ACCESS_HPP
#define PRIVATE_ACCESS_HPP

#include <utility>

namespace private_access
{
  template <typename C, typename T>
  struct member_wrapper
  {
    using type = T(C::*);
  };

  template <class T>
  struct proxy
  {
    static typename T::type value;
  };

  template <class T>
  typename T::type proxy<T>::value;

  template <class T, typename T::type AccessPointer>
  class make_proxy
  {
    struct Setter
    {
      Setter() { proxy<T>::value = AccessPointer; }
    };
    static Setter instance;
  };

  template <class T, typename T::type AccessPointer>
  typename make_proxy<T, AccessPointer>::Setter make_proxy<T, AccessPointer>::instance;

  template <typename Sig, class Instance>
  decltype(auto) member(Instance &&instance)
  {
    return instance.*(proxy<Sig>::value);
  }
}

#define PRIVATE_ACCESS_MEMBER_EX(accessor_name, cls, member, ret_type)   \
  struct accessor_name : ::private_access::member_wrapper<cls, ret_type> \
  {                                                                      \
  };                                                                     \
  template class ::private_access::make_proxy<accessor_name, &cls::member>;

#define PRIVATE_ACCESS_MEMBER(cls, member, ret_type) \
  PRIVATE_ACCESS_MEMBER_EX(cls##_##member, cls, member, ret_type)

#endif
