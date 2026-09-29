# Reverse proxy image of the nightly web deploy (dev/deploy/hetzner/): Caddy with
# the rate_limit handler of github.com/mholt/caddy-ratelimit compiled in, because
# stock Caddy has no rate limiting. dev/deploy/build-web-image builds it with the
# same tag as the app image, so a deploy and its rollback swap both together.
FROM caddy:2.11.4-builder@sha256:369218c81ca6d6af249981221b3a5c764d886dd5b058f51d144066de13f2418d AS builder

RUN xcaddy build v2.11.4 \
  --output /usr/bin/caddy \
  --with github.com/mholt/caddy-ratelimit@5625512f24f6f59d6f64fb3aafe5eecff0b286db

FROM caddy:2.11.4@sha256:14a9c00d4e833ebc2b65d36515b37bde3b73f0b323a2663aaafc88953d8c4e3f

COPY --from=builder /usr/bin/caddy /usr/bin/caddy
