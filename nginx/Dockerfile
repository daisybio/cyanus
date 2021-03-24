FROM debian:jessie
MAINTAINER mlist@health.sdu.dk
#set versions
ENV NGINX_VERSION 1.9.4

#install build dependencies
RUN apt-get update && apt-get install -y libpcre3-dev build-essential libssl-dev libgeoip-dev git wget

#add nginx user
RUN adduser --system --no-create-home --disabled-login --disabled-password --group nginx

#download and extract source
RUN cd /opt &&\
    git clone https://bitbucket.org/nginx-goodies/nginx-sticky-module-ng.git &&\
    wget http://nginx.org/download/nginx-$NGINX_VERSION.tar.gz -P/opt &&\
    tar -xvzf nginx-$NGINX_VERSION.tar.gz &&\
cd nginx-$NGINX_VERSION &&\
./configure --prefix=/opt/nginx \
	--user=nginx \
	--group=nginx \
	--with-http_ssl_module \
	--with-ipv6 \
        --with-http_geoip_module \
        --with-http_spdy_module \
        --add-module=/opt/nginx-sticky-module-ng \
	--error-log-path=/var/log/nginx/error.log \
	--http-log-path=/var/log/nginx/access.log &&\
make && make install

#clean up
RUN apt-get remove --purge -y build-essential && apt-get autoremove --purge -y && rm -rf /var/lib/apt/lists/*
EXPOSE 80 443

# forward request and error logs to docker log collector
RUN ln -sf /dev/stdout /var/log/nginx/access.log
RUN ln -sf /dev/stderr /var/log/nginx/error.log

#start nginx
CMD ["/opt/nginx/sbin/nginx", "-g", "daemon off;"]
