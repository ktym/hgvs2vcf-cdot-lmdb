# frozen_string_literal: true

source "https://rubygems.org"

gem "lmdb"              # embedded, file-persisted, mmap KV (text + coordinate index)
gem "sinatra"          # thin HTTP layer; swap for Roda if you prefer
gem "puma"             # concurrent app server
gem "rackup"           # Rack 3 CLI

group :development, :test do
  gem "minitest"
  gem "rack-test"
end
