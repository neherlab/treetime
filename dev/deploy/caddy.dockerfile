# Reverse proxy image of the nightly web deploy (dev/deploy/hetzner/): Caddy with
# the rate_limit handler of github.com/mholt/caddy-ratelimit compiled in, because
# stock Caddy has no rate limiting. dev/deploy/build-web-image builds it with the
# same tag as the app image, so a deploy and its rollback swap both together.
FROM caddy:2.11.7-builder@sha256:8910fe662f448b860e8d16ef1c9a4f8946d3193bdc3e971bf0c55c9a2cf246f6 AS builder

RUN xcaddy build v2.11.7 \
  --output /usr/bin/caddy \
  --with github.com/mholt/caddy-ratelimit@5625512f24f6f59d6f64fb3aafe5eecff0b286db

FROM caddy:2.11.7@sha256:685ae3d7c34d9b5f8caeddc5d14051db7bbd88ff1f5a2fa0c16bd8ca328fb7af

# The base image ships an older CA bundle.
RUN apk add --no-cache --upgrade "ca-certificates>=20261020-r0" "ca-certificates-bundle>=20261020-r0"

COPY --from=builder /usr/bin/caddy /usr/bin/caddy
